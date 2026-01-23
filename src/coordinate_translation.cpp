#include "pangenome_index/r-index.hpp"
#include "pangenome_index/sampled_tag_array.hpp"
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/simple_sds.hpp>
#include <gbwt/gbwt.h>
#include <gbwt/fast_locate.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gbz.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <sstream>
#include <limits>
#include <map>
#include <chrono>
#include <iomanip>

using namespace std;
using namespace std::chrono;

// Use explicit namespace for panindexer types to avoid conflict with gbwt::FastLocate
using panindexer::FastLocate;
using panindexer::SampledTagArray;

static bool debug = false;

// Structure to store tag information with source offsets
struct TagInfo {
    uint64_t tag_code;                    // Encoded (node_id, is_rev)
    vector<size_t> source_offsets;        // Offsets on source haplotype where this tag appears
    vector<size_t> source_bwt_positions; // BWT positions where source visits this tag
    vector<size_t> source_packed_positions; // Packed text positions for source visits
};

// Structure to store a visit to a node by a haplotype
struct NodeVisit {
    size_t seq_id;
    size_t offset;              // Offset on the haplotype
    size_t bwt_pos;             // BWT position
    size_t packed_pos;          // Packed text position (SA value)
    uint64_t tag_code;          // Tag code for this node
};

// Structure for coordinate translation result
struct TranslationResult {
    size_t source_offset;       // Offset on source haplotype
    size_t target_offset;       // Offset on target haplotype
    size_t target_seq_id;       // Target sequence ID (for fragments)
    uint64_t tag_code;          // Node ID at this position
};

// Structure to store fragment information
struct FragmentInfo {
    size_t seq_id;              // Sequence ID of this fragment
    size_t start_offset;        // Start offset on this fragment
    size_t end_offset;          // End offset on this fragment
    vector<size_t> anchor_offsets; // Anchor offsets where this fragment aligns
};

// Parse interval string like "1000..2000"
static pair<size_t, size_t> parse_interval(const string& interval_str) {
    size_t dot_pos = interval_str.find("..");
    if (dot_pos == string::npos) {
        throw runtime_error("Invalid interval format: expected 'start..end'");
    }
    size_t start = stoull(interval_str.substr(0, dot_pos));
    size_t end = stoull(interval_str.substr(dot_pos + 2));
    if (start > end) {
        throw runtime_error("Invalid interval: start > end");
    }
    return {start, end};
}

// Find sequence ID by direct path name (for generic paths like "x", "y", "chr1")
// Returns numeric_limits<size_t>::max() if not found
static size_t find_sequence_id_by_path_name(const gbwt::GBWT& gbwt_index,
                                             const string& path_name_str,
                                             bool debug_output = false) {
    if (!gbwt_index.hasMetadata()) {
        cerr << "Error: GBWT index does not have metadata" << endl;
        return numeric_limits<size_t>::max();
    }
    
    const gbwt::Metadata& metadata = gbwt_index.metadata;
    
    if (debug_output) {
        cerr << "Searching for path by name: '" << path_name_str << "'" << endl;
        cerr << "GBWT has " << metadata.paths() << " paths" << endl;
    }
    
    // For generic paths, the path name is typically stored as the contig name
    // First try to find a contig that matches the path name
    size_t contig_id = metadata.contig(path_name_str);
    
    if (contig_id < metadata.contigs()) {
        // Found a matching contig, now find paths with this contig
        if (debug_output) {
            cerr << "Found contig_id=" << contig_id << " for path name '" << path_name_str << "'" << endl;
        }
        
        vector<size_t> matching_paths;
        for (size_t path_id = 0; path_id < metadata.paths(); ++path_id) {
            gbwt::PathName pn = metadata.path(path_id);
            if (pn.contig == contig_id) {
                matching_paths.push_back(path_id);
                if (debug_output) {
                    cerr << "Found matching path: path_id=" << path_id 
                         << " (sample=" << pn.sample
                         << ", contig=" << pn.contig
                         << ", phase=" << pn.phase
                         << ", count=" << pn.count << ")" << endl;
                }
            }
        }
        
        if (!matching_paths.empty()) {
            if (matching_paths.size() > 1) {
                cerr << "Warning: Multiple paths (" << matching_paths.size() 
                     << ") found for contig '" << path_name_str << "'. Using first match (path_id=" 
                     << matching_paths[0] << ")" << endl;
            }
            
            if (debug_output) {
                cerr << "Resolved path name '" << path_name_str << "' to sequence ID: " << matching_paths[0] << endl;
            }
            return matching_paths[0];
        }
    }
    
    // If contig lookup failed, list available paths
    cerr << "Error: Path '" << path_name_str << "' not found" << endl;
    cerr << "Available contigs:" << endl;
    for (size_t i = 0; i < metadata.contigs(); ++i) {
        cerr << "  " << metadata.contig(i) << endl;
    }
    
    return numeric_limits<size_t>::max();
}

// Find sequence ID from path metadata (sample, contig, haplotype)
// Returns numeric_limits<size_t>::max() if not found
static size_t find_sequence_id_from_metadata(const gbwt::GBWT& gbwt_index,
                                              const string& sample_name,
                                              const string& contig_name,
                                              size_t haplotype,
                                              bool debug_output = false) {
    if (!gbwt_index.hasMetadata()) {
        cerr << "Error: GBWT index does not have metadata" << endl;
        return numeric_limits<size_t>::max();
    }
    
    const gbwt::Metadata& metadata = gbwt_index.metadata;
    
    if (debug_output) {
        cerr << "Searching for path: sample='" << sample_name 
             << "', contig='" << contig_name 
             << "', haplotype=" << haplotype << endl;
        cerr << "GBWT has " << metadata.paths() << " paths, "
             << metadata.samples() << " samples, "
             << metadata.contigs() << " contigs" << endl;
    }
    
    // Find sample ID from sample name
    size_t sample_id = metadata.sample(sample_name);
    if (sample_id >= metadata.samples()) {
        cerr << "Error: Sample '" << sample_name << "' not found in metadata" << endl;
        cerr << "Available samples:" << endl;
        for (size_t i = 0; i < metadata.samples(); ++i) {
            cerr << "  " << i << ": " << metadata.sample(i) << endl;
        }
        return numeric_limits<size_t>::max();
    }
    
    // Find contig ID from contig name
    size_t contig_id = metadata.contig(contig_name);
    if (contig_id >= metadata.contigs()) {
        cerr << "Error: Contig '" << contig_name << "' not found in metadata" << endl;
        cerr << "Available contigs:" << endl;
        for (size_t i = 0; i < metadata.contigs(); ++i) {
            cerr << "  " << i << ": " << metadata.contig(i) << endl;
        }
        return numeric_limits<size_t>::max();
    }
    
    if (debug_output) {
        cerr << "Found sample_id=" << sample_id << " for '" << sample_name << "'" << endl;
        cerr << "Found contig_id=" << contig_id << " for '" << contig_name << "'" << endl;
    }
    
    // Search through all paths to find matching one
    vector<size_t> matching_paths;
    for (size_t path_id = 0; path_id < metadata.paths(); ++path_id) {
        gbwt::PathName path_name = metadata.path(path_id);
        
        if (path_name.sample == sample_id && 
            path_name.contig == contig_id &&
            path_name.phase == haplotype) {
            matching_paths.push_back(path_id);
            if (debug_output) {
                cerr << "Found matching path: path_id=" << path_id 
                     << " (sample=" << path_name.sample
                     << ", contig=" << path_name.contig
                     << ", phase=" << path_name.phase
                     << ", count=" << path_name.count << ")" << endl;
            }
        }
    }
    
    if (matching_paths.empty()) {
        cerr << "Error: No path found matching sample='" << sample_name 
             << "', contig='" << contig_name 
             << "', haplotype=" << haplotype << endl;
        return numeric_limits<size_t>::max();
    }
    
    if (matching_paths.size() > 1) {
        cerr << "Warning: Multiple paths (" << matching_paths.size() 
             << ") found matching the criteria. Using first match (path_id=" 
             << matching_paths[0] << ")" << endl;
    }
    
    // The path_id in GBWT metadata corresponds to the sequence ID
    size_t seq_id = matching_paths[0];
    
    if (debug_output) {
        cerr << "Resolved to sequence ID: " << seq_id << endl;
    }
    
    return seq_id;
}

// Structure to hold path specification options for source or target
struct PathSpec {
    size_t seq_id = numeric_limits<size_t>::max();
    string path_name;
    string sample_name;
    string contig_name;
    size_t haplotype = 0;
    bool has_path_name = false;
    bool has_sample = false;
    bool has_contig = false;
    bool has_haplotype = false;
    bool reverse_strand = false;  // If true, use reverse complement strand (seq 2i+1 instead of 2i)
    
    bool use_path_name() const { return has_path_name; }
    bool use_metadata() const { return has_sample || has_contig; }
    bool use_seq_id() const { return seq_id != numeric_limits<size_t>::max(); }
    
    // Resolve the sequence ID using the provided GBWT index
    // For RLBWT, sequences are stored in pairs: seq 2i = forward, seq 2i+1 = reverse complement
    // When resolving from path name or metadata, we get the GBWT path_id and convert to RLBWT seq_id
    // Returns true on success, false on error
    bool resolve(const gbwt::GBWT& gbwt_index, const string& label, bool debug_output) {
        if (use_seq_id()) {
            // Already have seq_id (direct RLBWT seq_id), nothing to do
            // Note: when using --*-id directly, user provides the RLBWT seq_id
            return true;
        }
        
        size_t gbwt_path_id = numeric_limits<size_t>::max();
        
        if (use_path_name()) {
            gbwt_path_id = find_sequence_id_by_path_name(gbwt_index, path_name, debug_output);
            if (gbwt_path_id == numeric_limits<size_t>::max()) {
                cerr << "Error: Could not find " << label << " path '" << path_name << "'" << endl;
                return false;
            }
        } else if (use_metadata()) {
            gbwt_path_id = find_sequence_id_from_metadata(gbwt_index, sample_name, contig_name, haplotype, debug_output);
            if (gbwt_path_id == numeric_limits<size_t>::max()) {
                cerr << "Error: Could not resolve " << label << " path from metadata" << endl;
                return false;
            }
        } else {
            return false;
        }
        
        // Convert GBWT path_id to RLBWT seq_id
        // RLBWT stores both strands: seq 2i = forward, seq 2i+1 = reverse complement
        seq_id = 2 * gbwt_path_id + (reverse_strand ? 1 : 0);
        
        if (debug_output) {
            cerr << "Resolved " << label << ": GBWT path_id=" << gbwt_path_id 
                 << " -> RLBWT seq_id=" << seq_id 
                 << " (strand: " << (reverse_strand ? "reverse" : "forward") << ")" << endl;
            if (use_path_name()) {
                cerr << "  (from path name '" << path_name << "')" << endl;
            } else if (use_metadata()) {
                cerr << "  (from sample='" << sample_name 
                     << "', contig='" << contig_name 
                     << "', haplotype=" << haplotype << ")" << endl;
            }
        }
        
        return true;
    }
    
    // Validate that only one method is used
    bool validate(const string& label) const {
        int methods = (use_seq_id() ? 1 : 0) + (use_path_name() ? 1 : 0) + (use_metadata() ? 1 : 0);
        if (methods > 1) {
            cerr << "Error: Cannot mix --" << label << "-id, --" << label << "-path-name, and --" 
                 << label << "-sample/--" << label << "-contig options" << endl;
            return false;
        }
        if (use_metadata()) {
            if (!has_sample) {
                cerr << "Error: --" << label << "-sample is required when using --" << label << "-contig" << endl;
                return false;
            }
            if (!has_contig) {
                cerr << "Error: --" << label << "-contig is required when using --" << label << "-sample" << endl;
                return false;
            }
        }
        return true;
    }
    
    // Check if any specification method was provided
    bool is_specified() const {
        return use_seq_id() || use_path_name() || use_metadata();
    }
};

static void usage(const char* prog) {
    cerr << "Usage: " << prog << " <rlbwt_rindex.ri> <gbwt_rindex.ri> <sampled.tags> [options]" << endl;
    cerr << endl;
    cerr << "Required files:" << endl;
    cerr << "  rlbwt_rindex.ri         RLBWT R-index file" << endl;
    cerr << "  gbwt_rindex.ri          GBWT FastLocate R-index file (.ri)" << endl;
    cerr << "  sampled.tags            Sampled tag array file" << endl;
    cerr << endl;
    cerr << "Source specification (choose one):" << endl;
    cerr << "  --source-id ID                            Direct RLBWT sequence ID (already accounts for strand)" << endl;
    cerr << "  --source-path-name NAME [--source-reverse]" << endl;
    cerr << "                                            Specify by path name (e.g., 'x', 'chr1')" << endl;
    cerr << "  --source-sample NAME --source-contig NAME [--source-haplotype N] [--source-reverse]" << endl;
    cerr << "                                            Specify by metadata" << endl;
    cerr << endl;
    cerr << "Target specification (choose one):" << endl;
    cerr << "  --target-id ID                            Direct RLBWT sequence ID (already accounts for strand)" << endl;
    cerr << "  --target-path-name NAME [--target-reverse]" << endl;
    cerr << "                                            Specify by path name (e.g., 'x', 'chr1')" << endl;
    cerr << "  --target-sample NAME --target-contig NAME [--target-haplotype N] [--target-reverse]" << endl;
    cerr << "                                            Specify by metadata" << endl;
    cerr << endl;
    cerr << "Note on RLBWT sequence IDs:" << endl;
    cerr << "  RLBWT stores both strands: seq 2i = forward, seq 2i+1 = reverse complement" << endl;
    cerr << "  When using --*-path-name or --*-sample/--*-contig, the GBWT path_id is converted:" << endl;
    cerr << "    RLBWT seq_id = 2 * GBWT_path_id (forward) or 2 * GBWT_path_id + 1 (reverse)" << endl;
    cerr << "  Use --source-reverse or --target-reverse to select reverse complement strand" << endl;
    cerr << endl;
    cerr << "Required options:" << endl;
    cerr << "  --interval START..END   Interval on source haplotype (base offsets)" << endl;
    cerr << "  --gbwt-index FILE       GBWT index file (.gbwt) for LF mapping" << endl;
    cerr << "  --gbz FILE              GBZ file (.gbz) for GBWTGraph (node lengths)" << endl;
    cerr << endl;
    cerr << "Other options:" << endl;
    cerr << "  --debug                 Enable debug output" << endl;
    cerr << "  --no-debug              Disable debug output" << endl;
    cerr << endl;
    cerr << "Examples:" << endl;
    cerr << "  # Using direct RLBWT sequence IDs" << endl;
    cerr << "  " << prog << " idx.ri gbwt.ri tags.stags --source-id 0 --target-id 2 --interval 1000..2000 --gbwt-index graph.gbwt --gbz graph.gbz" << endl;
    cerr << endl;
    cerr << "  # Using path names (forward strand)" << endl;
    cerr << "  " << prog << " idx.ri gbwt.ri tags.stags --source-path-name x --target-path-name y --interval 1000..2000 --gbwt-index graph.gbwt --gbz graph.gbz" << endl;
    cerr << endl;
    cerr << "  # Using path names with reverse complement strand for target" << endl;
    cerr << "  " << prog << " idx.ri gbwt.ri tags.stags --source-path-name x --target-path-name y --target-reverse --interval 1000..2000 --gbwt-index graph.gbwt --gbz graph.gbz" << endl;
    cerr << endl;
    cerr << "  # Using metadata" << endl;
    cerr << "  " << prog << " idx.ri gbwt.ri tags.stags --source-sample HG002 --source-contig chr1 --source-haplotype 1 \\" << endl;
    cerr << "       --target-sample GRCh38 --target-contig chr1 --interval 1000..2000 --gbwt-index graph.gbwt --gbz graph.gbz" << endl;
}

// Find all tags (nodes) visited by the source haplotype in the given interval
vector<TagInfo> find_tags_in_interval(FastLocate& r_index, SampledTagArray& sampled,
                                       size_t source_seq_id, size_t seq_start, size_t seq_end) {
    vector<TagInfo> tags;
    unordered_map<uint64_t, TagInfo> tag_map;
    
    // Convert sequence interval to packed text positions
    size_t text_pos_i = r_index.pack(source_seq_id, seq_start);
    size_t text_pos_j = r_index.pack(source_seq_id, seq_end);
    
    if (debug) {
        cerr << "Finding tags in interval [" << seq_start << ", " << seq_end 
             << "] (text positions [" << text_pos_i << ", " << text_pos_j << "])" << endl;
    }
    
    // Find successor position on last at or after j
    r_index.ensure_last_rank();
    r_index.ensure_last_select();
    
    auto successor_result = r_index.last_successor(text_pos_j);
    size_t text_pos_x = successor_result.first;
    size_t rank_x = successor_result.second;
    
    size_t run_id = (rank_x < r_index.last_to_run.size()) 
        ? r_index.last_to_run[rank_x] : 0;
    size_t bwt_pos_x = r_index.bwt_end_position_of_run(run_id);
    
    if (debug) {
        cerr << "Found marked text position x=" << text_pos_x << " (rank=" << rank_x 
             << ", run_id=" << run_id << ", BWT_pos=" << bwt_pos_x << ")" << endl;
    }
    
    // Check if text_pos_x is after the sequence end
    size_t bwt_seq_end = source_seq_id;  // BWT position of $ at end of sequence source_seq_id
    size_t text_pos_seq_end;
    {
        size_t rindex_run_id = 0;
        size_t run_start_pos = 0;
        r_index.run_id_and_offset_at(bwt_seq_end, rindex_run_id, run_start_pos);
        text_pos_seq_end = r_index.getSample(rindex_run_id);
        // Navigate from run start to bwt_seq_end
        for (size_t p = run_start_pos; p < bwt_seq_end; ++p) {
            text_pos_seq_end = r_index.locateNext(text_pos_seq_end);
        }
    }
    
    if (debug) {
        cerr << "Sequence end check: text_pos_x=" << text_pos_x 
             << ", text_pos_seq_end=" << text_pos_seq_end << " (BWT[" << bwt_seq_end << "])" << endl;
    }
    
    // Check if text_pos_x (BWT rank from last array) is greater than sequence end
    if (text_pos_x > text_pos_seq_end) {
        cerr << "Error: BWT rank calculated from last array (text_pos_x=" << text_pos_x 
             << ") is greater than sequence end (text_pos_seq_end=" << text_pos_seq_end << ")" << endl;
        return tags;  // Return empty tags vector
    }
    
    // Iterate backwards through the interval
    sampled.ensure_run_rank();
    sampled.ensure_run_select();
    
    size_t current_text_pos = text_pos_x;
    size_t current_bwt_pos = bwt_pos_x;
    
    while (current_text_pos >= text_pos_i) {
        if (current_text_pos <= text_pos_j) {
            // Get tag for current BWT position
            size_t sampled_run_id = sampled.run_id_at(current_bwt_pos);
            uint64_t tag_val = sampled.run_value(sampled_run_id);

            if (debug) {
                cerr << "  Text pos " << current_text_pos << " -> BWT pos " << current_bwt_pos 
                     << " -> tag=" << tag_val << endl;
            }
            
            if (tag_val != 0) {
                // Unpack to get offset
                size_t seq_id_found = r_index.seqId(current_text_pos);
                size_t offset_from_start = r_index.seqOffset(current_text_pos);
                
                // Only process positions in our source sequence and interval
                if (seq_id_found == source_seq_id && 
                    offset_from_start >= seq_start && offset_from_start <= seq_end) {
                    
                    if (tag_map.find(tag_val) == tag_map.end()) {
                        tag_map[tag_val] = TagInfo{tag_val, {}, {}, {}};
                    }
                    
                    tag_map[tag_val].source_offsets.push_back(offset_from_start);
                    tag_map[tag_val].source_bwt_positions.push_back(current_bwt_pos);
                    tag_map[tag_val].source_packed_positions.push_back(current_text_pos);
                }
            }
        }
        
        // Move to previous text position
        if (current_text_pos > text_pos_i) {
            current_text_pos--;
            current_bwt_pos = r_index.LF(current_bwt_pos);
        } else {
            break;
        }
    }
    
    // Convert map to vector and sort by first source offset (earliest in interval)
    for (auto& kv : tag_map) {
        tags.push_back(kv.second);
    }
    
    // Sort by first source offset (prefer nodes near start of interval)
    sort(tags.begin(), tags.end(), [](const TagInfo& a, const TagInfo& b) {
        if (a.source_offsets.empty() || b.source_offsets.empty()) return false;
        return a.source_offsets[0] < b.source_offsets[0];
    });
    
    if (debug) {
        cerr << "Found " << tags.size() << " unique tags in source interval" << endl;
    }
    // print all the tags
    if (debug) {
        cerr << "Tags found in source interval:" << endl;
        for (const auto& tag : tags) {
            cerr << "Tag code: " << tag.tag_code << endl;
            cerr << "Source offsets: [";
            for (size_t i = 0; i < tag.source_offsets.size(); i++) {
                if (i > 0) cerr << ", ";
                cerr << tag.source_offsets[i];
            }
            cerr << "]" << endl;
        }
    }
    
    return tags;
}

// Find all sequences (and their offsets) that pass through a given tag
vector<NodeVisit> find_sequences_for_tag(FastLocate& r_index, SampledTagArray& sampled,
                                          uint64_t tag_code) {
    vector<NodeVisit> visits;
    
    if (debug) {
        cerr << "[find_sequences_for_tag] Searching for tag_code=" << tag_code << endl;
    }
    
    const auto& wm = sampled.values();
    size_t total_occ = wm.rank(wm.size(), tag_code);
    
    if (debug) {
        cerr << "  wm.size()=" << wm.size() << ", total_occ=" << total_occ << endl;
    }
    
    if (total_occ == 0) {
        if (debug) {
            cerr << "  No occurrences found, returning empty" << endl;
        }
        return visits;
    }
    
    bool first_run_is_gap = sampled.is_first_run_gap();
    if (debug) {
        cerr << "  first_run_is_gap=" << first_run_is_gap << endl;
    }
    
    for (size_t j = 1; j <= total_occ; ++j) {
        if (debug) {
            cerr << "  Processing occurrence " << j << "/" << total_occ << endl;
        }
        
        size_t tag_run_id = wm.select(j, tag_code);
        if (debug) {
            cerr << "    tag_run_id=" << tag_run_id << endl;
        }
        
        size_t bwt_run_id = 2 * tag_run_id + 1 - static_cast<size_t>(first_run_is_gap);
        if (debug) {
            cerr << "    bwt_run_id=" << bwt_run_id << endl;
        }
        
        auto run_span = sampled.run_span(bwt_run_id);
        size_t locate_start = run_span.first;
        size_t locate_end = run_span.second;
        
        if (debug) {
            cerr << "    run_span=[" << locate_start << ", " << locate_end 
                 << "], size=" << (locate_end - locate_start + 1) << endl;
        }
        
        // IMPORTANT: Skip positions >= r_index.bwt_size() (those are endmarkers)
        // The sampled tag array includes endmarkers, but RLBWT r-index doesn't
        size_t r_index_size = r_index.bwt_size();
        if (locate_start >= r_index_size) {
            if (debug) {
                cerr << "    Skipping entire run: locate_start=" << locate_start 
                     << " >= r_index.bwt_size()=" << r_index_size << " (endmarker positions)" << endl;
            }
            continue;  // Skip this entire run
        }
        
        // Clamp locate_end to valid range
        if (locate_end >= r_index_size) {
            if (debug) {
                cerr << "    Clamping locate_end from " << locate_end << " to " << (r_index_size - 1) 
                     << " (r_index.bwt_size()=" << r_index_size << ")" << endl;
            }
            locate_end = r_index_size - 1;
        }
        
        // Get initial packed position
        size_t packed_pos;
        {
            size_t rindex_run_id = 0;
            size_t run_start_pos = 0;
            
            if (debug) {
                cerr << "    Getting initial packed position at locate_start=" << locate_start << endl;
                cerr << "      r_index.bwt_size()=" << r_index.bwt_size() << endl;
                cerr << "      r_index.size() (num runs)=" << r_index.size() << endl;
            }
            
            r_index.run_id_and_offset_at(locate_start, rindex_run_id, run_start_pos);
            
            if (debug) {
                cerr << "      rindex_run_id=" << rindex_run_id << ", run_start_pos=" << run_start_pos << endl;
            }
            
            packed_pos = r_index.getSample(rindex_run_id);
            
            if (debug) {
                cerr << "      initial packed_pos=" << packed_pos << endl;
            }
            
            // Navigate from run start to locate_start
            size_t nav_steps = locate_start - run_start_pos;
            if (debug) {
                cerr << "      Need to navigate " << nav_steps << " steps from run_start_pos=" 
                     << run_start_pos << " to locate_start=" << locate_start << endl;
            }
            
            for (size_t p = run_start_pos; p < locate_start; ++p) {
                if (debug) {
                    cerr << "      Before locateNext: p=" << p << ", packed_pos=" << packed_pos << endl;
                }
                packed_pos = r_index.locateNext(packed_pos);
                if (debug) {
                    cerr << "      After locateNext: p=" << p << ", packed_pos=" << packed_pos << endl;
                }
            }
            
            if (debug) {
                cerr << "      After navigation: packed_pos=" << packed_pos << endl;
            }
        }
        
        // Process all positions in the run
        for (size_t pos = locate_start; pos <= locate_end; ++pos) {
            if (debug) {
                cerr << "    Processing pos=" << pos << " (offset=" << (pos - locate_start) 
                     << "/" << (locate_end - locate_start) << ")" << endl;
            }
            if (debug && ((pos - locate_start) % 100 == 0 || pos == locate_start || pos == locate_end)) {
                cerr << "    Processing pos=" << pos << " (offset=" << (pos - locate_start) 
                     << "/" << (locate_end - locate_start) << ")" << endl;
            }
            
            if (pos > locate_start) {
                if (debug) {
                    cerr << "      Before locateNext: pos=" << pos << ", packed_pos=" << packed_pos << endl;
                    cerr << "        Checking if packed_pos is valid (< bwt_size=" << r_index.bwt_size() << ")..." << endl;
                }
                
                // if (packed_pos >= r_index.bwt_size()) {
                //     cerr << "ERROR: packed_pos=" << packed_pos << " is >= bwt_size=" << r_index.bwt_size() << endl;
                //     cerr << "  This should not happen! packed_pos should be < bwt_size" << endl;
                //     return visits;  // Return what we have so far
                // }
                
                packed_pos = r_index.locateNext(packed_pos);
                
                if (debug) {
                    cerr << "      After locateNext: pos=" << pos << ", packed_pos=" << packed_pos << endl;
                }
                if (debug && ((pos - locate_start) % 100 == 0 || pos == locate_end)) {
                    cerr << "      After locateNext: packed_pos=" << packed_pos << endl;
                }
            }
            
            // Unpack to get (seq_id, offset)
            if (debug && ((pos - locate_start) % 100 == 0 || pos == locate_start || pos == locate_end)) {
                cerr << "      About to unpack packed_pos=" << packed_pos << endl;
            }
            auto pr = r_index.unpack(packed_pos);

            size_t seq_id_found = pr.first;
            size_t offset_from_start = pr.second;
            
            if (debug && ((pos - locate_start) % 100 == 0 || pos == locate_start || pos == locate_end)) {
                cerr << "      Unpacked: seq_id=" << seq_id_found << ", offset=" << offset_from_start << endl;
            }
            
            NodeVisit visit;
            visit.seq_id = seq_id_found;
            visit.offset = offset_from_start;
            visit.bwt_pos = pos;
            visit.packed_pos = packed_pos;
            visit.tag_code = tag_code;
            
            visits.push_back(visit);
            // std::cerr << "done finding sequences for tag" << endl;
        }
        
        if (debug) {
            cerr << "    Finished occurrence " << j << ", total visits so far: " << visits.size() << endl;
        }
    }
    
    if (debug) {
        cerr << "  Total visits found: " << visits.size() << endl;
    }
    
    return visits;
}

// Find all common nodes between source and target haplotypes (handles fragments)
// Returns vector of anchor points: (source_offset, target_offset, tag_code, source_packed_pos, target_packed_pos, target_seq_id)
vector<tuple<size_t, size_t, uint64_t, size_t, size_t, size_t>> find_all_common_nodes(
    FastLocate& r_index, SampledTagArray& sampled,
    const vector<TagInfo>& source_tags, size_t target_seq_id) {
    
    vector<tuple<size_t, size_t, uint64_t, size_t, size_t, size_t>> anchors;
    
    for (const auto& tag_info : source_tags) {
        // Find all sequences that pass through this tag
        vector<NodeVisit> all_visits = find_sequences_for_tag(r_index, sampled, tag_info.tag_code);
        
        // Find visits by target haplotype (may include multiple fragments)
        vector<NodeVisit> target_visits;
        for (const auto& visit : all_visits) {
            if (visit.seq_id == target_seq_id) {
                target_visits.push_back(visit);
            }
        }
        
        if (!target_visits.empty()) {
            // Found common node(s)! Match visits by SA value ordering
            
            for (size_t src_idx = 0; src_idx < tag_info.source_offsets.size(); src_idx++) {
                size_t source_offset = tag_info.source_offsets[src_idx];
                size_t source_packed_pos = tag_info.source_packed_positions[src_idx];
                
                // Find target visit with closest SA value
                size_t best_target_idx = 0;
                size_t min_diff = numeric_limits<size_t>::max();
                
                for (size_t tgt_idx = 0; tgt_idx < target_visits.size(); tgt_idx++) {
                    size_t diff = (source_packed_pos > target_visits[tgt_idx].packed_pos) 
                        ? (source_packed_pos - target_visits[tgt_idx].packed_pos)
                        : (target_visits[tgt_idx].packed_pos - source_packed_pos);
                    
                    if (diff < min_diff) {
                        min_diff = diff;
                        best_target_idx = tgt_idx;
                    }
                }
                
                size_t target_offset = target_visits[best_target_idx].offset;
                size_t target_packed_pos = target_visits[best_target_idx].packed_pos;
                size_t target_fragment_seq_id = target_visits[best_target_idx].seq_id;
                
                if (debug) {
                    uint64_t decoded = tag_info.tag_code - 1;
                    int64_t node_id = (decoded >> 1) + 1;
                    bool is_rev = (decoded & 1) != 0;
                    cerr << "Found common node: node_id=" << node_id 
                         << ", is_rev=" << is_rev
                         << ", source_offset=" << source_offset
                         << ", target_seq_id=" << target_fragment_seq_id
                         << ", target_offset=" << target_offset << endl;
                }
                
                anchors.push_back(make_tuple(source_offset, target_offset, tag_info.tag_code, 
                                            source_packed_pos, target_packed_pos, target_fragment_seq_id));
            }
        }
    }
    
    return anchors;
}

// Decode tag code to get node ID and orientation
pair<int64_t, bool> decode_tag(uint64_t tag_code) {
    uint64_t decoded = tag_code - 1;
    int64_t node_id = (decoded >> 1) + 1;
    bool is_rev = (decoded & 1) != 0;
    return {node_id, is_rev};
}

// Convert node offset to base offset for a given GBWT sequence
// Returns the cumulative base offset at the start of the node at node_offset
// Requires GBWTGraph to get node lengths
size_t node_offset_to_base_offset(
    const gbwt::GBWT& gbwt_index,
    const gbwtgraph::GBWTGraph& graph,
    size_t seq_id,
    size_t node_offset) {
    
    // Extract the path for this sequence
    gbwt::vector_type path = gbwt_index.extract(gbwt::Path::encode(seq_id, false));
    
    if (path.empty()) {
        cerr << "Warning: Sequence " << seq_id << " not found in GBWT" << endl;
        return 0;
    }
    
    size_t cumulative_bases = 0;
    size_t current_node_offset = 0;
    
    for (gbwt::node_type node : path) {
        if (node == gbwt::ENDMARKER) break;
        if (current_node_offset >= node_offset) break;
        
        // Get node handle from graph
        gbwtgraph::handle_t handle = graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node));
        size_t node_length = graph.get_length(handle);
        
        cumulative_bases += node_length;
        current_node_offset++;
    }
    
    return cumulative_bases;
}

// Convert base offset to node offset for a given GBWT sequence
// Returns the node offset (number of nodes from start) that contains the base offset
// Requires GBWTGraph to get node lengths
size_t base_offset_to_node_offset(
    const gbwt::GBWT& gbwt_index,
    const gbwtgraph::GBWTGraph& graph,
    size_t seq_id,
    size_t base_offset) {
    
    if (debug) {
        cerr << "[base_offset_to_node_offset] seq_id=" << seq_id 
             << ", base_offset=" << base_offset << endl;
    }
    
    // Extract the path for this sequence
    gbwt::vector_type path = gbwt_index.extract(gbwt::Path::encode(seq_id, false));
    
    if (path.empty()) {
        cerr << "Warning: Sequence " << seq_id << " not found in GBWT" << endl;
        return 0;
    }
    
    if (debug) {
        cerr << "  Path length: " << path.size() << " nodes" << endl;
    }
    
    size_t cumulative_bases = 0;
    size_t node_offset = 0;
    
    for (gbwt::node_type node : path) {
        if (node == gbwt::ENDMARKER) break;
        
        // Get node handle from graph
        gbwtgraph::handle_t handle = graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node));
        size_t node_length = graph.get_length(handle);
        
        if (debug) {
            cerr << "  node_offset=" << node_offset << ", node_id=" << gbwt::Node::id(node)
                 << ", cumulative_bases=[" << cumulative_bases << ", " 
                 << (cumulative_bases + node_length) << ")" << endl;
        }
        
        // Check if base_offset is within this node
        // cumulative_bases is the start of current node, cumulative_bases + node_length is the end
        if (base_offset < cumulative_bases + node_length) {
            // Base offset is within this node
            if (debug) {
                cerr << "  Found! base_offset " << base_offset << " is in node_offset " << node_offset << endl;
            }
            return node_offset;
        }
        
        cumulative_bases += node_length;
        node_offset++;
    }
    
    // Base offset is beyond the path, return last node offset
    if (debug) {
        cerr << "  Warning: base_offset " << base_offset << " is beyond path end (cumulative_bases=" 
             << cumulative_bases << "), returning last node_offset=" << (node_offset > 0 ? node_offset - 1 : 0) << endl;
    }
    return node_offset > 0 ? node_offset - 1 : 0;
}

// Structure to hold first and last common nodes
struct CommonNodes {
    size_t first_source_offset;      // GBWT node offset
    size_t first_target_offset;      // GBWT node offset
    size_t first_source_base;        // RLBWT base offset
    size_t first_target_base;        // RLBWT base offset
    uint64_t first_tag_code;
    size_t last_source_offset;       // GBWT node offset
    size_t last_target_offset;      // GBWT node offset
    size_t last_source_base;         // RLBWT base offset
    size_t last_target_base;        // RLBWT base offset
    uint64_t last_tag_code;
    bool found;
};

// Structure to store a node visit along a path during traversal
struct PathNode {
    uint64_t tag_code;               // Tag code for this node
    gbwt::node_type node;            // GBWT node
    size_t node_offset;              // Node offset from start of path
    size_t base_offset;               // Base offset from start of path
    gbwt::edge_type gbwt_position;   // GBWT edge position for LF mapping
};

// Helper function to check if a tag is a common node and extract offsets
// use_largest_offset controls offset selection:
//   false: smallest RLBWT base + largest GBWT offset (both = earliest) -> for FIRST node
//   true:  largest RLBWT base + smallest GBWT offset (both = latest) -> for LAST node
bool check_common_node(
    const gbwt::FastLocate& gbwt_fast_locate,
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const TagInfo& tag_info,
    size_t source_seq_id, size_t target_seq_id,
    size_t& source_offset, size_t& target_offset,
    size_t& source_base, size_t& target_base,
    bool use_largest_offset = false) {
    
    // Decode tag to get GBWT node
    auto [node_id, is_rev] = decode_tag(tag_info.tag_code);
    gbwt::node_type node = gbwt::Node::encode(node_id, is_rev);
    
    if (debug) {
        cerr << "  Checking tag_code=" << tag_info.tag_code 
             << " (node_id=" << node_id << ", is_rev=" << is_rev << ")" << endl;
    }
    
    // Use GBWT FastLocate to find all paths visiting this node
    std::vector<gbwt::size_type> sa_values = gbwt_fast_locate.decompressSA(node);
    
    if (debug) {
        cerr << "    Found " << sa_values.size() << " path occurrences on this node" << endl;
    }
    
    // Find source and target visits
    // Note: In GBWT, larger offset = earlier in path (opposite of RLBWT)
    vector<pair<gbwt::size_type, size_t>> source_visits;  // (sa_value, path_offset)
    vector<pair<gbwt::size_type, size_t>> target_visits;  // (sa_value, path_offset)
    
    for (size_t i = 0; i < sa_values.size(); i++) {
        gbwt::size_type seq_id = gbwt_fast_locate.seqId(sa_values[i]);
        gbwt::size_type seq_offset = gbwt_fast_locate.seqOffset(sa_values[i]);

        if (debug) {
            cerr << "    SA value: " << sa_values[i] << " (seq_id=" << seq_id << ", seq_offset=" << seq_offset << ")" << endl;
        }
        
        if (seq_id == source_seq_id) {
            source_visits.push_back({sa_values[i], seq_offset});
        }
        if (seq_id == target_seq_id) {
            target_visits.push_back({sa_values[i], seq_offset});
        }
    }
    
    if (debug) {
        cerr << "    Source visits: " << source_visits.size() << endl;
        cerr << "    Target visits: " << target_visits.size() << endl;
    }
    
    if (source_visits.empty() || target_visits.empty()) {
        return false;
    }
    
    // Found common node! Match visits based on offsets
    // Sort visits by offset (larger offset = earlier position in path)
    sort(source_visits.begin(), source_visits.end(),
         [](const pair<gbwt::size_type, size_t>& a, const pair<gbwt::size_type, size_t>& b) {
             return a.second > b.second;  // Larger offset first (earlier in path)
         });
    sort(target_visits.begin(), target_visits.end(),
         [](const pair<gbwt::size_type, size_t>& a, const pair<gbwt::size_type, size_t>& b) {
             return a.second > b.second;  // Larger offset first (earlier in path)
         });
    
    // Select GBWT offset based on whether we want first or last common node
    // IMPORTANT: In GBWT, LARGER offset = EARLIER in path (opposite of RLBWT!)
    if (use_largest_offset) {
        // use_largest_offset = true: LAST common node
        // Select SMALLEST GBWT offset (latest in path)
        source_offset = source_visits[source_visits.size() - 1].second;  // last index = smallest offset
        target_offset = target_visits[target_visits.size() - 1].second;
        if (debug) {
            cerr << "    Selected SMALLEST GBWT offsets (latest in path, for last common node): source=" 
                 << source_offset << ", target=" << target_offset << endl;
        }
    } else {
        // use_largest_offset = false: FIRST common node
        // Select LARGEST GBWT offset (earliest in path)
        source_offset = source_visits[0].second;  // index 0 = largest offset after sorting
        target_offset = target_visits[0].second;
        if (debug) {
            cerr << "    Selected LARGEST GBWT offsets (earliest in path, for first common node): source=" 
                 << source_offset << ", target=" << target_offset << endl;
        }
    }
    
    // Get base offsets from RLBWT r-index for this tag
    // Collect all source and target visits from RLBWT
    vector<NodeVisit> source_visits_rlbwt = find_sequences_for_tag(rlbwt_rindex, sampled, tag_info.tag_code);
    
    // Separate source and target visits
    vector<size_t> source_base_offsets;
    vector<size_t> target_base_offsets;
    
    for (const auto& visit : source_visits_rlbwt) {
        if (visit.seq_id == source_seq_id) {
            source_base_offsets.push_back(visit.offset);
        } else if (visit.seq_id == target_seq_id) {
            target_base_offsets.push_back(visit.offset);
        }
    }
    
    if (source_base_offsets.empty() || target_base_offsets.empty()) {
        return false;  // Not common to both
    }
    
    // Sort offsets (ascending order)
    sort(source_base_offsets.begin(), source_base_offsets.end());
    sort(target_base_offsets.begin(), target_base_offsets.end());
    
    // For RLBWT: smaller offset = earlier in sequence, larger offset = later in sequence
    if (use_largest_offset) {
        // use_largest_offset = true: LAST common node
        // Select largest RLBWT base offset (latest in sequence)
        source_base = source_base_offsets.back();
        target_base = target_base_offsets.back();
    } else {
        // use_largest_offset = false: FIRST common node
        // Select smallest RLBWT base offset (earliest in sequence)
        source_base = source_base_offsets[0];
        target_base = target_base_offsets[0];
    }

    if (debug) {
        cerr << "    RLBWT base offsets:" << endl;
        if (use_largest_offset) {
            cerr << "      Selected LARGEST RLBWT offsets (latest in sequence, for last common node):" << endl;
        } else {
            cerr << "      Selected SMALLEST RLBWT offsets (earliest in sequence, for first common node):" << endl;
        }
        cerr << "      Source base: " << source_base << endl;
        cerr << "      Target base: " << target_base << endl;
    }

    
    return true;
}

// Find the first (smallest tag_code) and last (largest tag_code) common nodes
// between source and target haplotypes using GBWT FastLocate
// Optimized: finds first by iterating from beginning, finds last by iterating from end
// Also looks up base offsets from RLBWT r-index
// Returns CommonNodes structure with both first and last common nodes
CommonNodes find_first_and_last_common_nodes_gbwt(
    const gbwt::FastLocate& gbwt_fast_locate,
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const vector<TagInfo>& source_tags, 
    size_t source_seq_id, size_t target_seq_id) {
    
    CommonNodes result;
    result.found = false;
    
    // Sort tags by their position in the source path (by first base offset, smallest to largest)
    // This ensures "first" means earliest in the path and "last" means latest in the path
    vector<TagInfo> sorted_tags = source_tags;
    sort(sorted_tags.begin(), sorted_tags.end(), 
         [](const TagInfo& a, const TagInfo& b) {
             // Sort by tag_code (ascending order)
             return a.tag_code < b.tag_code;
         });
    
    if (debug) {
        cerr << "Finding first and last common nodes (optimized search):" << endl;
        cerr << "  Total tags: " << sorted_tags.size() << endl;
        cerr << "  Tags sorted by source base offset:" << endl;
        for (size_t i = 0; i < sorted_tags.size(); i++) {
            auto [node_id, is_rev] = decode_tag(sorted_tags[i].tag_code);
            size_t offset = sorted_tags[i].source_offsets.empty() ? 0 : sorted_tags[i].source_offsets[0];
            cerr << "    [" << i << "] tag_code=" << sorted_tags[i].tag_code 
                 << " (node_id=" << node_id << "), base_offset=" << offset << endl;
        }
    }
    
    // Step 1: Find first common node by iterating from beginning (earliest in path)
    if (debug) {
        cerr << "  Searching for FIRST common node (earliest in path)..." << endl;
    }
    
    for (size_t i = 0; i < sorted_tags.size(); i++) {
        const auto& tag_info = sorted_tags[i];
        size_t source_offset, target_offset, source_base, target_base;
        
        // For first common node: smallest RLBWT offset + largest GBWT offset (both = earliest)
        if (check_common_node(gbwt_fast_locate, rlbwt_rindex, sampled, tag_info,
                              source_seq_id, target_seq_id,
                              source_offset, target_offset, source_base, target_base,
                              false)) {  // false = smallest RLBWT + largest GBWT (earliest)
            // Found first common node!
            result.first_source_offset = source_offset;
            result.first_target_offset = target_offset;
            result.first_source_base = source_base;
            result.first_target_base = target_base;
            result.first_tag_code = tag_info.tag_code;
            result.found = true;
            
            if (debug) {
                cerr << "  ✓ Found FIRST common node (earliest in path)!" << endl;
                cerr << "    (smallest RLBWT base, largest GBWT offset)" << endl;
                cerr << "    Source: GBWT_offset=" << source_offset 
                     << ", RLBWT_base=" << source_base << endl;
                cerr << "    Target: GBWT_offset=" << target_offset
                     << ", RLBWT_base=" << target_base << endl;
            }
            break;
        }
    }
    
    if (!result.found) {
        if (debug) {
            cerr << "  No common nodes found" << endl;
        }
        return result;
    }
    
    // Step 2: Find last common node by iterating from end (latest in path)
    if (debug) {
        cerr << "  Searching for LAST common node (latest in path)..." << endl;
    }
    
    for (size_t i = sorted_tags.size(); i > 0; i--) {
        const auto& tag_info = sorted_tags[i - 1];
        size_t source_offset, target_offset, source_base, target_base;
        
        // For last common node: largest RLBWT offset + smallest GBWT offset (both = latest)
        if (check_common_node(gbwt_fast_locate, rlbwt_rindex, sampled, tag_info,
                              source_seq_id, target_seq_id,
                              source_offset, target_offset, source_base, target_base,
                              true)) {  // true = largest RLBWT + smallest GBWT (latest)
            // Found last common node!
            result.last_source_offset = source_offset;
            result.last_target_offset = target_offset;
            result.last_source_base = source_base;
            result.last_target_base = target_base;
            result.last_tag_code = tag_info.tag_code;
            
            if (debug) {
                cerr << "  ✓ Found LAST common node (latest in path)!" << endl;
                cerr << "    (largest RLBWT base, smallest GBWT offset)" << endl;
                cerr << "    Source: GBWT_offset=" << source_offset 
                     << ", RLBWT_base=" << source_base << endl;
                cerr << "    Target: GBWT_offset=" << target_offset
                     << ", RLBWT_base=" << target_base << endl;
            }
            break;
        }
    }
    
    if (debug) {
        cerr << "\n  Summary:" << endl;
        cerr << "    First common node: tag_code=" << result.first_tag_code << endl;
        cerr << "      Source: GBWT_node_offset=" << result.first_source_offset 
             << ", RLBWT_base_offset=" << result.first_source_base << endl;
        cerr << "      Target: GBWT_node_offset=" << result.first_target_offset 
             << ", RLBWT_base_offset=" << result.first_target_base << endl;
        cerr << "    Last common node: tag_code=" << result.last_tag_code << endl;
        cerr << "      Source: GBWT_node_offset=" << result.last_source_offset 
             << ", RLBWT_base_offset=" << result.last_source_base << endl;
        cerr << "      Target: GBWT_node_offset=" << result.last_target_offset 
             << ", RLBWT_base_offset=" << result.last_target_base << endl;
    }
    
    return result;
}

// Trace paths using GBWT LF mapping (one node at a time)
// Uses GBWT's efficient forward traversal to map positions between source and target haplotypes
// FastLocate is used to find initial positions, then GBWT index LF() is used for traversal
unordered_map<size_t, vector<size_t>> build_offset_mapping_with_gbwt_lf(
    const gbwt::GBWT& gbwt_index,
    const gbwt::FastLocate& gbwt_fast_locate,
    const gbwtgraph::GBWTGraph& graph,
    size_t source_seq_id, size_t source_start, size_t source_end,
    size_t target_seq_id, size_t anchor_source_offset, size_t anchor_target_offset,
    size_t anchor_source_base, size_t anchor_target_base,
    uint64_t anchor_tag_code, size_t last_common_source_offset, uint64_t last_common_tag_code) {
    
    unordered_map<size_t, vector<size_t>> offset_map;
    
    if (debug) {
        cerr << "Building offset mapping using GBWT LF mapping:" << endl;
        cerr << "  Anchor source node offset: " << anchor_source_offset << endl;
        cerr << "  Anchor target node offset: " << anchor_target_offset << endl;
    }
    
    // Decode anchor node and encode as GBWT node
    auto [anchor_node_id, anchor_is_rev] = decode_tag(anchor_tag_code);
    gbwt::node_type anchor_node = gbwt::Node::encode(anchor_node_id, anchor_is_rev);
    
    if (debug) {
        cerr << "  Anchor node: id=" << anchor_node_id 
             << ", is_rev=" << anchor_is_rev 
             << ", encoded=" << anchor_node << endl;
    }
    
    // Decompress SA for just the anchor node (efficient - only this node's occurrences)
    std::vector<gbwt::size_type> sa_values = gbwt_fast_locate.decompressSA(anchor_node);
    
    if (debug) {
        cerr << "  SA values for anchor node: " << sa_values.size() << " occurrences" << endl;
    }
    
    // Find the offsets (indices) for source and target sequences in the SA result
    // The index in sa_values IS the offset needed for LF mapping
    gbwt::size_type source_offset_in_node = gbwt::invalid_offset();
    gbwt::size_type target_offset_in_node = gbwt::invalid_offset();
    
    for (size_t i = 0; i < sa_values.size(); i++) {
        gbwt::size_type seq_id = gbwt_fast_locate.seqId(sa_values[i]);
        gbwt::size_type seq_offset = gbwt_fast_locate.seqOffset(sa_values[i]);
        
        if (seq_id == source_seq_id && seq_offset == anchor_source_offset) {
            source_offset_in_node = i;
            if (debug) {
                cerr << "  Found source at index " << i << " in node's SA" << endl;
            }
        }
        if (seq_id == target_seq_id && seq_offset == anchor_target_offset) {
            target_offset_in_node = i;
            if (debug) {
                cerr << "  Found target at index " << i << " in node's SA" << endl;
            }
        }
    }
    
    // TODO: Handle case where multiple paths from same sequence pass through anchor node
    // For now, we assume only one occurrence per sequence at this node
    if (source_offset_in_node == gbwt::invalid_offset()) {
        cerr << "Warning: Source sequence " << source_seq_id 
             << " not found at anchor node (offset " << anchor_source_offset << ")" << endl;
        return offset_map;
    }
    if (target_offset_in_node == gbwt::invalid_offset()) {
        cerr << "Warning: Target sequence " << target_seq_id 
             << " not found at anchor node (offset " << anchor_target_offset << ")" << endl;
        return offset_map;
    }
    
    // Create GBWT positions for LF mapping
    // edge_type = (node, offset_within_node's_record)
    gbwt::edge_type source_pos = {anchor_node, source_offset_in_node};
    gbwt::edge_type target_pos = {anchor_node, target_offset_in_node};
    
    if (debug) {
        cerr << "  Starting positions for LF traversal:" << endl;
        cerr << "    Source: (node=" << source_pos.first << ", offset=" << source_pos.second << ")" << endl;
        cerr << "    Target: (node=" << target_pos.first << ", offset=" << target_pos.second << ")" << endl;
    }
    
    // Track current offsets in each path (in terms of node count from start of path)
    size_t source_path_offset = anchor_source_offset;
    size_t target_path_offset = anchor_target_offset;
    
    // Track cumulative base offsets (starting from RLBWT base offsets, then adding node lengths)
    size_t source_base_offset = anchor_source_base;
    size_t target_base_offset = anchor_target_base;
    
    // Start with anchor point (using base offsets from RLBWT)
    if (anchor_source_base >= source_start && anchor_source_base <= source_end) {
        offset_map[anchor_source_base].push_back(anchor_target_base);
    }
    
    if (debug) {
        cerr << "  Starting from anchor (base offsets from RLBWT):" << endl;
        cerr << "    Anchor: source_node=" << anchor_source_offset << ", source_base=" << anchor_source_base << endl;
        cerr << "    Anchor: target_node=" << anchor_target_offset << ", target_base=" << anchor_target_base << endl;
        cerr << "    Input interval: [" << source_start << ", " << source_end << "] bases" << endl;
    }
    
    // ========== FORWARD TRAVERSAL using GBWT LF mapping ==========
    // GBWT's LF mapping moves forward along the path (one node at a time)
    // Unlike standard FM-index LF which moves backward, GBWT's LF moves forward because:
    //   1. GBWT stores paths as sequences of nodes (not bases)
    //   2. Each node's record stores outgoing edges (next nodes)
    //   3. LF(pos) returns the next node position along the path
    //
    // IMPORTANT: Paths can diverge and converge again
    // Example: Source: 1→2→3, Target: 1→4→5→6→2→3
    // We traverse each path INDEPENDENTLY, store all nodes, then map at common nodes
    
    // Store paths as vectors of PathNode entries
    vector<PathNode> source_path;
    vector<PathNode> target_path;
    
    // Add anchor node to both paths
    PathNode anchor_src_node;
    anchor_src_node.tag_code = anchor_tag_code;
    anchor_src_node.node = anchor_node;
    anchor_src_node.node_offset = anchor_source_offset;
    anchor_src_node.base_offset = anchor_source_base;
    anchor_src_node.gbwt_position = source_pos;
    source_path.push_back(anchor_src_node);
    
    PathNode anchor_tgt_node;
    anchor_tgt_node.tag_code = anchor_tag_code;
    anchor_tgt_node.node = anchor_node;
    anchor_tgt_node.node_offset = anchor_target_offset;
    anchor_tgt_node.base_offset = anchor_target_base;
    anchor_tgt_node.gbwt_position = target_pos;
    target_path.push_back(anchor_tgt_node);
    
    if (debug) {
        cerr << "  Starting independent path traversal..." << endl;
        cerr << "  Anchor node: tag_code=" << anchor_tag_code << endl;
        cerr << "    Source: node_offset=" << anchor_source_offset << ", base_offset=" << anchor_source_base << endl;
        cerr << "    Target: node_offset=" << anchor_target_offset << ", base_offset=" << anchor_target_base << endl;
    }
    
    // ========== Traverse source path independently ==========
    gbwt::edge_type src_pos = source_pos;
    size_t src_node_offset = anchor_source_offset;
    size_t src_base_offset = anchor_source_base;
    
    if (debug) {
        cerr << "  Traversing source path..." << endl;
    }
    
    while (true) {
        gbwt::edge_type next_src_pos = gbwt_index.LF(src_pos);
        
        if (next_src_pos.first == gbwt::ENDMARKER) {
            if (debug) {
                cerr << "    Source path stopped: reached ENDMARKER after " << source_path.size() << " nodes" << endl;
            }
            break;
        }
        
        // Update base offset by adding length of node we just left
        gbwtgraph::handle_t curr_src_handle = graph.get_handle(
            gbwt::Node::id(src_pos.first), gbwt::Node::is_reverse(src_pos.first));
        size_t node_len = graph.get_length(curr_src_handle);
        src_base_offset += node_len;
        src_node_offset++;
        
        if (debug) {
            cerr << "    Moved forward: added node_len=" << node_len 
                 << ", new src_base_offset=" << src_base_offset << endl;
        }
        
        // Get tag code for next node
        gbwt::node_type next_src_node = next_src_pos.first;
        int64_t next_src_node_id = gbwt::Node::id(next_src_node);
        bool next_src_is_rev = gbwt::Node::is_reverse(next_src_node);
        uint64_t src_decoded = ((static_cast<uint64_t>(next_src_node_id) - 1) << 1) | (next_src_is_rev ? 1 : 0);
        uint64_t next_src_tag_code = src_decoded + 1;
        
        // Check stop conditions
        if (next_src_tag_code > last_common_tag_code) {
            if (debug) {
                cerr << "    Source path stopped: tag_code=" << next_src_tag_code 
                     << " > last_common_tag_code=" << last_common_tag_code 
                     << " after " << source_path.size() << " nodes" << endl;
            }
            break;
        }
        
        // Stop if we've gone past the source interval
        // base_offset is the START of the node, so stop if start > source_end
        if (src_base_offset > source_end) {
            if (debug) {
                cerr << "    Source path stopped: src_base_offset=" << src_base_offset 
                     << " > source_end=" << source_end 
                     << " after " << source_path.size() << " nodes" << endl;
            }
            break;
        }
        
        // Store this node (only if within reasonable bounds)
        PathNode src_node;
        src_node.tag_code = next_src_tag_code;
        src_node.node = next_src_node;
        src_node.node_offset = src_node_offset;
        src_node.base_offset = src_base_offset;
        src_node.gbwt_position = next_src_pos;
        source_path.push_back(src_node);
        
        if (debug) {
            cerr << "    Added source node: tag_code=" << next_src_tag_code 
                 << ", node_offset=" << src_node_offset 
                 << ", base_offset=" << src_base_offset << endl;
        }
        
        src_pos = next_src_pos;
    }
    
    // ========== Traverse target path independently ==========
    gbwt::edge_type tgt_pos = target_pos;
    size_t tgt_node_offset = anchor_target_offset;
    size_t tgt_base_offset = anchor_target_base;
    
    if (debug) {
        cerr << "  Traversing target path..." << endl;
    }
    
    while (true) {
        gbwt::edge_type next_tgt_pos = gbwt_index.LF(tgt_pos);
        
        if (next_tgt_pos.first == gbwt::ENDMARKER) {
            if (debug) {
                cerr << "    Target path stopped: reached ENDMARKER after " << target_path.size() << " nodes" << endl;
            }
            break;
        }
        
        // Update base offset by adding length of node we just left
        gbwtgraph::handle_t curr_tgt_handle = graph.get_handle(
            gbwt::Node::id(tgt_pos.first), gbwt::Node::is_reverse(tgt_pos.first));
        size_t node_len = graph.get_length(curr_tgt_handle);
        tgt_base_offset += node_len;
        tgt_node_offset++;
        
        if (debug) {
            cerr << "    Moved forward: added node_len=" << node_len 
                 << ", new tgt_base_offset=" << tgt_base_offset << endl;
        }
        
        // Get tag code for next node
        gbwt::node_type next_tgt_node = next_tgt_pos.first;
        int64_t next_tgt_node_id = gbwt::Node::id(next_tgt_node);
        bool next_tgt_is_rev = gbwt::Node::is_reverse(next_tgt_node);
        uint64_t tgt_decoded = ((static_cast<uint64_t>(next_tgt_node_id) - 1) << 1) | (next_tgt_is_rev ? 1 : 0);
        uint64_t next_tgt_tag_code = tgt_decoded + 1;
        
        // Check stop condition
        if (next_tgt_tag_code > last_common_tag_code) {
            if (debug) {
                cerr << "    Target path stopped: tag_code=" << next_tgt_tag_code 
                     << " > last_common_tag_code=" << last_common_tag_code 
                     << " after " << target_path.size() << " nodes" << endl;
            }
            break;
        }
        
        // Store this node
        PathNode tgt_node;
        tgt_node.tag_code = next_tgt_tag_code;
        tgt_node.node = next_tgt_node;
        tgt_node.node_offset = tgt_node_offset;
        tgt_node.base_offset = tgt_base_offset;
        tgt_node.gbwt_position = next_tgt_pos;
        target_path.push_back(tgt_node);
        
        if (debug) {
            cerr << "    Added target node: tag_code=" << next_tgt_tag_code 
                 << ", node_offset=" << tgt_node_offset 
                 << ", base_offset=" << tgt_base_offset << endl;
        }
        
        tgt_pos = next_tgt_pos;
    }
    
    if (debug) {
        cerr << "  Path traversal completed:" << endl;
        cerr << "    Source path: " << source_path.size() << " nodes" << endl;
        cerr << "    Target path: " << target_path.size() << " nodes" << endl;
    }
    
    // ========== Find common nodes and map positions ==========
    // Create a map from tag_code to target path nodes for efficient lookup
    unordered_map<uint64_t, vector<size_t>> target_nodes_by_tag;
    for (size_t i = 0; i < target_path.size(); i++) {
        target_nodes_by_tag[target_path[i].tag_code].push_back(i);
    }
    
    if (debug) {
        cerr << "  Finding common nodes and mapping positions..." << endl;
    }
    
    // Iterate through source path and find matching nodes in target path
    for (const auto& src_node : source_path) {
        // Only map positions within the source interval
        if (src_node.base_offset < source_start || src_node.base_offset > source_end) {
            continue;
        }
        
        // Check if this tag_code exists in target path
        auto it = target_nodes_by_tag.find(src_node.tag_code);
        if (it != target_nodes_by_tag.end()) {
            // Found common node!
            // NOTE: If there are multiple occurrences in target, we map to ALL of them
            // This handles cases where paths diverge and reconverge
            
            if (it->second.size() > 1 && debug) {
                cerr << "    Common node with MULTIPLE target occurrences: tag_code=" << src_node.tag_code 
                     << " (" << it->second.size() << " occurrences)" << endl;
            }
            
            for (size_t tgt_idx : it->second) {
                const PathNode& tgt_node = target_path[tgt_idx];
                
                // Store ALL mappings (multiple target offsets per source offset)
                offset_map[src_node.base_offset].push_back(tgt_node.base_offset);
                
                if (debug) {
                    cerr << "    Common node: tag_code=" << src_node.tag_code << endl;
                    cerr << "      Source: base_offset=" << src_node.base_offset << ", node_offset=" << src_node.node_offset << endl;
                    cerr << "      Target: base_offset=" << tgt_node.base_offset << ", node_offset=" << tgt_node.node_offset << endl;
                    cerr << "      -> Mapped: " << src_node.base_offset << " -> " << tgt_node.base_offset << endl;
                }
            }
        }
    }
    
    if (debug) {
        cerr << "  Mapping complete: " << offset_map.size() << " positions mapped" << endl;
    }
    
    return offset_map;
}

// Legacy wrapper using panindexer::FastLocate (for compatibility)
// TODO: Remove this once all callers are updated to use GBWT directly
unordered_map<size_t, vector<size_t>> build_offset_mapping_with_gbwt_rindex(
    FastLocate& rlbwt_rindex, FastLocate& gbwt_rindex,
    size_t source_seq_id, size_t source_start, size_t source_end,
    size_t target_seq_id, size_t anchor_source_offset, size_t anchor_target_offset,
    uint64_t anchor_tag_code) {
    
    cerr << "Warning: Using legacy build_offset_mapping_with_gbwt_rindex. "
         << "Please update to use build_offset_mapping_with_gbwt_lf with GBWT objects." << endl;
    
    // Return just the anchor point mapping as fallback
    unordered_map<size_t, vector<size_t>> offset_map;
    offset_map[anchor_source_offset].push_back(anchor_target_offset);
    return offset_map;
}

// Trace coordinates forward from anchor point using GBWT LF mapping
vector<TranslationResult> trace_coordinates_gbwt(
    const gbwt::GBWT& gbwt_index,
    const gbwt::FastLocate& gbwt_fast_locate,
    const gbwtgraph::GBWTGraph& graph,
    size_t source_seq_id, size_t source_start, size_t source_end,
    size_t target_seq_id, size_t anchor_source_offset, size_t anchor_target_offset,
    size_t anchor_source_base, size_t anchor_target_base,
    uint64_t anchor_tag_code, size_t last_common_source_offset, uint64_t last_common_tag_code) {
    
    vector<TranslationResult> translations;
    
    if (debug) {
        cerr << "Tracing coordinates from anchor using GBWT LF mapping:" << endl;
        cerr << "  Anchor source offset: " << anchor_source_offset << endl;
        cerr << "  Anchor target offset: " << anchor_target_offset << endl;
    }
    
    // Build offset mapping using GBWT LF-based tracing
    unordered_map<size_t, vector<size_t>> offset_map = build_offset_mapping_with_gbwt_lf(
        gbwt_index, gbwt_fast_locate, graph, source_seq_id, source_start, source_end,
        target_seq_id, anchor_source_offset, anchor_target_offset,
        anchor_source_base, anchor_target_base,
        anchor_tag_code, last_common_source_offset, last_common_tag_code);
    
    if (debug) {
        cerr << "\n=== Final Offset Mapping ===" << endl;
        cerr << "Total source positions mapped: " << offset_map.size() << endl;
        for (const auto& [src, tgt_vec] : offset_map) {
            cerr << "  " << src << " -> [ ";
            for (size_t i = 0; i < tgt_vec.size(); i++) {
                if (i > 0) cerr << ", ";
                cerr << tgt_vec[i];
            }
            cerr << " ]" << endl;
        }
        cerr << "========================\n" << endl;
    }
    
    // Convert to translation results (create one result per source->target mapping)
    for (size_t src_off = source_start; src_off <= source_end; src_off++) {
        if (offset_map.find(src_off) != offset_map.end()) {
            // For each source offset, create a result for EACH target offset
            const vector<size_t>& target_offsets = offset_map[src_off];
            for (size_t tgt_off : target_offsets) {
                TranslationResult result;
                result.source_offset = src_off;
                result.target_seq_id = target_seq_id;
                result.target_offset = tgt_off;
                result.tag_code = anchor_tag_code;
                translations.push_back(result);
            }
        } else {
            // No mapping found for this source offset
            TranslationResult result;
            result.source_offset = src_off;
            result.target_seq_id = target_seq_id;
            result.target_offset = 0;
            result.tag_code = anchor_tag_code;
            translations.push_back(result);
            
            if (debug) {
                cerr << "Warning: Source offset " << src_off << " not mapped to target" << endl;
            }
        }
    }
    
    return translations;
}

// Trace coordinates forward and backward from anchor point using legacy r-index
vector<TranslationResult> trace_coordinates(
    FastLocate& rlbwt_rindex, FastLocate& gbwt_rindex,
    size_t source_seq_id, size_t source_start, size_t source_end,
    size_t target_seq_id, size_t anchor_source_offset, size_t anchor_target_offset,
    uint64_t anchor_tag_code) {
    
    vector<TranslationResult> translations;
    
    if (debug) {
        cerr << "Tracing coordinates from anchor using legacy r-index:" << endl;
        cerr << "  Anchor source offset: " << anchor_source_offset << endl;
        cerr << "  Anchor target offset: " << anchor_target_offset << endl;
    }
    
    // Build offset mapping using legacy r-index tracing
    unordered_map<size_t, vector<size_t>> offset_map = build_offset_mapping_with_gbwt_rindex(
        rlbwt_rindex, gbwt_rindex, source_seq_id, source_start, source_end,
        target_seq_id, anchor_source_offset, anchor_target_offset,
        anchor_tag_code);
    
    // Convert to translation results (create one result per source->target mapping)
    for (size_t src_off = source_start; src_off <= source_end; src_off++) {
        if (offset_map.find(src_off) != offset_map.end()) {
            // For each source offset, create a result for EACH target offset
            const vector<size_t>& target_offsets = offset_map[src_off];
            for (size_t tgt_off : target_offsets) {
                TranslationResult result;
                result.source_offset = src_off;
                result.target_seq_id = target_seq_id;
                result.target_offset = tgt_off;
                result.tag_code = anchor_tag_code;
                translations.push_back(result);
            }
        } else {
            // No mapping found for this source offset
            TranslationResult result;
            result.source_offset = src_off;
            result.target_seq_id = target_seq_id;
            result.target_offset = 0;
            result.tag_code = anchor_tag_code;
            translations.push_back(result);
        }
    }
    
    return translations;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        usage(argv[0]);
        return 1;
    }
    
    string rlbwt_rindex_file = argv[1];
    string gbwt_rindex_file = argv[2];
    string sampled_tags_file = argv[3];
    
    // Parse arguments
    PathSpec source_spec;
    PathSpec target_spec;
    pair<size_t, size_t> interval = {0, 0};
    bool has_interval = false;
    string gbwt_index_file = "";
    string gbz_file = "";
    
    for (int i = 4; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--interval" && i + 1 < argc) {
            interval = parse_interval(argv[++i]);
            has_interval = true;
        } 
        // Source specification options
        else if (arg == "--source-id" && i + 1 < argc) {
            source_spec.seq_id = stoull(argv[++i]);
        } else if (arg == "--source-path-name" && i + 1 < argc) {
            source_spec.path_name = argv[++i];
            source_spec.has_path_name = true;
        } else if (arg == "--source-sample" && i + 1 < argc) {
            source_spec.sample_name = argv[++i];
            source_spec.has_sample = true;
        } else if (arg == "--source-contig" && i + 1 < argc) {
            source_spec.contig_name = argv[++i];
            source_spec.has_contig = true;
        } else if (arg == "--source-haplotype" && i + 1 < argc) {
            source_spec.haplotype = stoull(argv[++i]);
            source_spec.has_haplotype = true;
        }
        // Target specification options
        else if (arg == "--target-id" && i + 1 < argc) {
            target_spec.seq_id = stoull(argv[++i]);
        } else if (arg == "--target-path-name" && i + 1 < argc) {
            target_spec.path_name = argv[++i];
            target_spec.has_path_name = true;
        } else if (arg == "--target-sample" && i + 1 < argc) {
            target_spec.sample_name = argv[++i];
            target_spec.has_sample = true;
        } else if (arg == "--target-contig" && i + 1 < argc) {
            target_spec.contig_name = argv[++i];
            target_spec.has_contig = true;
        } else if (arg == "--target-haplotype" && i + 1 < argc) {
            target_spec.haplotype = stoull(argv[++i]);
            target_spec.has_haplotype = true;
        } else if (arg == "--target-reverse") {
            target_spec.reverse_strand = true;
        }
        // Source strand option
        else if (arg == "--source-reverse") {
            source_spec.reverse_strand = true;
        }
        // Other options
        else if (arg == "--gbwt-index" && i + 1 < argc) {
            gbwt_index_file = argv[++i];
        } else if (arg == "--gbz" && i + 1 < argc) {
            gbz_file = argv[++i];
        } else if (arg == "--debug") {
            debug = true;
        } else if (arg == "--no-debug") {
            debug = false;
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return 0;
        } else {
            cerr << "Unknown argument: " << arg << endl;
            usage(argv[0]);
            return 1;
        }
    }
    
    // Validate source specification
    if (!source_spec.validate("source")) {
        usage(argv[0]);
        return 1;
    }
    if (!source_spec.is_specified()) {
        cerr << "Error: Source haplotype must be specified (--source-id, --source-path-name, or --source-sample/--source-contig)" << endl;
        usage(argv[0]);
        return 1;
    }
    
    // Validate target specification
    if (!target_spec.validate("target")) {
        usage(argv[0]);
        return 1;
    }
    if (!target_spec.is_specified()) {
        cerr << "Error: Target haplotype must be specified (--target-id, --target-path-name, or --target-sample/--target-contig)" << endl;
        usage(argv[0]);
        return 1;
    }
    
    if (!has_interval) {
        cerr << "Error: --interval is required" << endl;
        usage(argv[0]);
        return 1;
    }
    
    if (gbwt_index_file.empty()) {
        cerr << "Error: --gbwt-index is required" << endl;
        usage(argv[0]);
        return 1;
    }
    
    if (gbz_file.empty()) {
        cerr << "Error: --gbz is required for GBWTGraph (needed for node lengths)" << endl;
        usage(argv[0]);
        return 1;
    }
    
    size_t seq_start = interval.first;
    size_t seq_end = interval.second;
    
    auto start_total = high_resolution_clock::now();
    
    // Load RLBWT r-index
    FastLocate rlbwt_rindex;
    {
        if (debug) cerr << "Loading RLBWT r-index: " << rlbwt_rindex_file << "..." << endl;
        ifstream rin(rlbwt_rindex_file, ios::binary);
        if (!rin) { 
            cerr << "Cannot open RLBWT r-index: " << rlbwt_rindex_file << endl; 
            return 1; 
        }
        rlbwt_rindex.load_encoded(rin);
    }
    
    // Load GBWT index - required for LF mapping
    // Note: FastLocate needs the GBWT index it was built from for LF operations
    if (debug) cerr << "Loading GBWT index: " << gbwt_index_file << "..." << endl;
    std::unique_ptr<gbwt::GBWT> gbwt_index_ptr;
    {
        ifstream gin_idx(gbwt_index_file, ios::binary);
        if (!gin_idx) {
            cerr << "Error: Cannot open GBWT index file: " << gbwt_index_file << endl;
            return 1;
        }
        try {
            gbwt_index_ptr = std::make_unique<gbwt::GBWT>();
            sdsl::load_from_file(*gbwt_index_ptr, gbwt_index_file);
            if (debug) cerr << "  GBWT index loaded. Sequences: " << gbwt_index_ptr->sequences() << endl;
        } catch (const std::exception& e) {
            cerr << "Error loading GBWT index: " << e.what() << endl;
            return 1;
        }
    }
    
    // Load GBWT FastLocate (r-index) from .ri file
    if (debug) cerr << "Loading GBWT r-index (FastLocate): " << gbwt_rindex_file << "..." << endl;
    std::unique_ptr<gbwt::FastLocate> gbwt_rindex_ptr;
    {
        ifstream gin(gbwt_rindex_file, ios::binary);
        if (!gin) { 
            cerr << "Cannot open GBWT r-index: " << gbwt_rindex_file << endl; 
            return 1; 
        }
        
        // Check file size
        gin.seekg(0, ios::end);
        size_t file_size = gin.tellg();
        gin.seekg(0, ios::beg);
        if (debug) cerr << "  File size: " << file_size << " bytes" << endl;
        
        if (file_size == 0) {
            cerr << "Error: GBWT r-index file is empty" << endl;
            return 1;
        }
        
        try {
            gbwt_rindex_ptr = std::make_unique<gbwt::FastLocate>();
            gbwt_rindex_ptr->load(gin);
            if (debug) cerr << "  GBWT r-index loaded successfully" << endl;
        } catch (const std::exception& e) {
            cerr << "Error loading GBWT r-index: " << e.what() << endl;
            return 1;
        } catch (...) {
            cerr << "Unknown error loading GBWT r-index (possibly segfault)" << endl;
            return 1;
        }
    }
    
    // Associate FastLocate with the GBWT index (required for decompressSA and other operations)
    // FastLocate stores a pointer to the GBWT index, which must be set after loading from file
    if (gbwt_index_ptr && gbwt_rindex_ptr) {
        gbwt_rindex_ptr->setGBWT(*gbwt_index_ptr);
        if (debug) cerr << "  FastLocate associated with GBWT index" << endl;
    }
    
    // Resolve source and target path specifications to sequence IDs
    // This requires the GBWT index to be loaded
    if (!source_spec.use_seq_id() || !target_spec.use_seq_id()) {
        if (!gbwt_index_ptr) {
            cerr << "Error: GBWT index is required for path name resolution" << endl;
            return 1;
        }
        
        if (!source_spec.resolve(*gbwt_index_ptr, "source", debug)) {
            return 1;
        }
        
        if (!target_spec.resolve(*gbwt_index_ptr, "target", debug)) {
            return 1;
        }
    }
    
    size_t source_seq_id = source_spec.seq_id;
    size_t target_seq_id = target_spec.seq_id;
    
    if (debug) {
        cerr << "Resolved sequence IDs:" << endl;
        cerr << "  Source: " << source_seq_id;
        if (source_spec.use_path_name()) {
            cerr << " (from path name '" << source_spec.path_name << "')";
        } else if (source_spec.use_metadata()) {
            cerr << " (from sample='" << source_spec.sample_name 
                 << "', contig='" << source_spec.contig_name 
                 << "', haplotype=" << source_spec.haplotype << ")";
        }
        cerr << endl;
        
        cerr << "  Target: " << target_seq_id;
        if (target_spec.use_path_name()) {
            cerr << " (from path name '" << target_spec.path_name << "')";
        } else if (target_spec.use_metadata()) {
            cerr << " (from sample='" << target_spec.sample_name 
                 << "', contig='" << target_spec.contig_name 
                 << "', haplotype=" << target_spec.haplotype << ")";
        }
        cerr << endl;
    }
    
    // Load GBWTGraph from GBZ file (required for node lengths to convert base/node offsets)
    if (debug) cerr << "Loading GBWTGraph from GBZ: " << gbz_file << "..." << endl;
    gbwtgraph::GBZ gbz;
    {
        try {
            sdsl::simple_sds::load_from(gbz, gbz_file);
            if (debug) cerr << "  GBWTGraph loaded successfully" << endl;
        } catch (const std::exception& e) {
            cerr << "Error loading GBZ file: " << e.what() << endl;
            return 1;
        }
    }
    const gbwtgraph::GBWTGraph& graph = gbz.graph;
    
    // Load sampled tag array
    SampledTagArray sampled;
    {
        if (debug) cerr << "Loading sampled tag array: " << sampled_tags_file << "..." << endl;
        ifstream sin(sampled_tags_file, ios::binary);
        if (!sin) { 
            cerr << "Cannot open sampled.tags: " << sampled_tags_file << endl; 
            return 1; 
        }
        sampled.load(sin);
    }
    
    if (debug) {
        cerr << "\n=== All files loaded successfully ===" << endl;
        cerr << "RLBWT R-index: loaded" << endl;
        if (gbwt_index_ptr && gbwt_rindex_ptr) {
            cerr << "GBWT index: loaded (" << gbwt_index_ptr->sequences() << " sequences)" << endl;
            cerr << "GBWT FastLocate: loaded" << endl;
        }
        cerr << "Sampled tag array: loaded" << endl;
        cerr << "=====================================\n" << endl;
    }
    
    // Step 1: Find all tags in source interval (using RLBWT r-index)
    auto start_find_tags = high_resolution_clock::now();
    vector<TagInfo> source_tags = find_tags_in_interval(rlbwt_rindex, sampled, source_seq_id, seq_start, seq_end);
    auto end_find_tags = high_resolution_clock::now();
    auto duration_find_tags = duration_cast<milliseconds>(end_find_tags - start_find_tags);
    
    if (source_tags.empty()) {
        cerr << "No tags found in source interval" << endl;
        return 1;
    }
    
    // Step 2: Find first and last common nodes with target haplotype using GBWT FastLocate
    // Sorts tags by tag_code (smallest to largest) and uses decompressSA() to find matching paths
    auto start_find_common = high_resolution_clock::now();
    CommonNodes common_nodes;
    
    if (gbwt_index_ptr && gbwt_rindex_ptr) {
        common_nodes = find_first_and_last_common_nodes_gbwt(*gbwt_rindex_ptr, rlbwt_rindex, sampled, source_tags, source_seq_id, target_seq_id);
    } else {
        cerr << "Error: GBWT index and FastLocate required for finding common nodes" << endl;
        return 1;
    }
    
    auto end_find_common = high_resolution_clock::now();
    auto duration_find_common = duration_cast<milliseconds>(end_find_common - start_find_common);
    
    if (!common_nodes.found) {
        cerr << "No common nodes found between source and target haplotypes" << endl;
        return 1;
    }
    
    size_t anchor_source_offset = common_nodes.first_source_offset;
    size_t anchor_target_offset = common_nodes.first_target_offset;
    size_t anchor_source_base = common_nodes.first_source_base;
    size_t anchor_target_base = common_nodes.first_target_base;
    uint64_t anchor_tag_code = common_nodes.first_tag_code;
    size_t last_common_source_offset = common_nodes.last_source_offset;
    uint64_t last_common_tag_code = common_nodes.last_tag_code;
    
    // Step 3: Trace coordinates using GBWT LF mapping
    auto start_trace = high_resolution_clock::now();
    vector<TranslationResult> translations;
    
    if (gbwt_index_ptr && gbwt_rindex_ptr) {
        // Use GBWT LF mapping for forward path traversal
        // FastLocate finds initial positions, GBWT index LF() does the traversal
        // Traversal stops when reaching the last common node
        translations = trace_coordinates_gbwt(
            *gbwt_index_ptr, *gbwt_rindex_ptr, graph, source_seq_id, seq_start, seq_end,
            target_seq_id, anchor_source_offset, anchor_target_offset,
            anchor_source_base, anchor_target_base, anchor_tag_code,
            last_common_source_offset, last_common_tag_code);
    } else {
        cerr << "Error: Both GBWT index and GBWT r-index are required" << endl;
        return 1;
    }
    
    auto end_trace = high_resolution_clock::now();
    auto duration_trace = duration_cast<milliseconds>(end_trace - start_trace);
    
    // Output results
    cout << "Coordinate Translation Results:" << endl;
    cout << "Source haplotype:" << endl;
    cout << "  Sequence ID: " << source_seq_id << endl;
    if (source_spec.use_path_name()) {
        cout << "  Path name: " << source_spec.path_name << endl;
        cout << "  Strand: " << (source_spec.reverse_strand ? "reverse" : "forward") << endl;
    } else if (source_spec.use_metadata()) {
        cout << "  Sample: " << source_spec.sample_name << endl;
        cout << "  Contig: " << source_spec.contig_name << endl;
        cout << "  Haplotype: " << source_spec.haplotype << endl;
        cout << "  Strand: " << (source_spec.reverse_strand ? "reverse" : "forward") << endl;
    }
    cout << "  Interval: " << seq_start << ".." << seq_end << endl;
    cout << "Target haplotype:" << endl;
    cout << "  Sequence ID: " << target_seq_id << endl;
    if (target_spec.use_path_name()) {
        cout << "  Path name: " << target_spec.path_name << endl;
        cout << "  Strand: " << (target_spec.reverse_strand ? "reverse" : "forward") << endl;
    } else if (target_spec.use_metadata()) {
        cout << "  Sample: " << target_spec.sample_name << endl;
        cout << "  Contig: " << target_spec.contig_name << endl;
        cout << "  Haplotype: " << target_spec.haplotype << endl;
        cout << "  Strand: " << (target_spec.reverse_strand ? "reverse" : "forward") << endl;
    }
    cout << "First common node (anchor): tag_code=" << anchor_tag_code << endl;
    cout << "  Source: node_offset=" << anchor_source_offset 
         << ", base_offset=" << anchor_source_base << endl;
    cout << "  Target: node_offset=" << anchor_target_offset 
         << ", base_offset=" << anchor_target_base << endl;
    cout << "Last common node (stop): tag_code=" << last_common_tag_code << endl;
    cout << "  Source: node_offset=" << last_common_source_offset 
         << ", base_offset=" << common_nodes.last_source_base << endl;
    cout << endl;
    
    // Check for fragments
    auto all_anchors = find_all_common_nodes(rlbwt_rindex, sampled, source_tags, target_seq_id);
    
    cout << "Translations:" << endl;
    cout << "Source_Offset\tTarget_Seq_ID\tTarget_Offset" << endl;
    for (const auto& trans : translations) {
        cout << trans.source_offset << "\t" << trans.target_seq_id << "\t" 
             << trans.target_offset << endl;
    }
    
    auto end_total = high_resolution_clock::now();
    auto duration_total = duration_cast<milliseconds>(end_total - start_total);
    
    if (debug) {
        cerr << "\n=== Timing Statistics ===" << endl;
        cerr << "  Find tags in interval: " << duration_find_tags.count() << " ms" << endl;
        cerr << "  Find common node: " << duration_find_common.count() << " ms" << endl;
        cerr << "  Trace coordinates: " << duration_trace.count() << " ms" << endl;
        cerr << "  Total time: " << duration_total.count() << " ms" << endl;
    }
    
    return 0;
}