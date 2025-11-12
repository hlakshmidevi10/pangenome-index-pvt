#include "pangenome_index/r-index.hpp"
#include "pangenome_index/sampled_tag_array.hpp"
#include <sdsl/wavelet_trees.hpp>
#include <gbwtgraph/utils.h>
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

using namespace std;
using namespace panindexer;
using namespace gbwtgraph;

struct TagResult {
    vector<size_t> query_offsets; // offsets within the query interval (relative to l)
    unordered_map<size_t, vector<size_t>> sequence_offsets; // seq_id -> list of offsets in that sequence
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

static void usage(const char* prog) {
    cerr << "Usage: " << prog << " <r_index.ri> <sampled.tags> <gbz_file> [--sample SAMPLE] [--contig CONTIG] [--haplotype H] [--interval START..END]" << endl;
    cerr << "  OR:   " << prog << " <r_index.ri> <sampled.tags> <gbz_file> --seq-id ID --interval START..END" << endl;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        usage(argv[0]);
        return 1;
    }
    
    string r_index_file = argv[1];
    string sampled_tags_file = argv[2];
    string gbz_file = argv[3];
    
    // Parse optional arguments
    string sample, contig;
    size_t haplotype = 0;
    size_t seq_id = numeric_limits<size_t>::max();
    pair<size_t, size_t> interval = {0, 0};
    bool has_interval = false;
    
    for (int i = 4; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--sample" && i + 1 < argc) {
            sample = argv[++i];
        } else if (arg == "--contig" && i + 1 < argc) {
            contig = argv[++i];
        } else if (arg == "--haplotype" && i + 1 < argc) {
            haplotype = stoull(argv[++i]);
        } else if (arg == "--interval" && i + 1 < argc) {
            interval = parse_interval(argv[++i]);
            has_interval = true;
        } else if (arg == "--seq-id" && i + 1 < argc) {
            seq_id = stoull(argv[++i]);
        } else {
            cerr << "Unknown argument: " << arg << endl;
            usage(argv[0]);
            return 1;
        }
    }
    
    if (!has_interval) {
        cerr << "Error: --interval is required" << endl;
        usage(argv[0]);
        return 1;
    }
    
    // Load GBZ file
    GBZ gbz;
    {
        cerr << "Loading GBZ file: " << gbz_file << "..." << endl;
        ifstream gbz_in(gbz_file, ios::binary);
        if (!gbz_in) {
            cerr << "Cannot open GBZ file: " << gbz_file << endl;
            return 1;
        }
        gbz_in.close();
    }
    sdsl::simple_sds::load_from(gbz, gbz_file);
    cerr << "GBZ file loaded successfully. Total sequences: " << gbz.index.sequences() << endl;
    
    // Determine sequence ID
    if (seq_id == numeric_limits<size_t>::max()) {
        if (sample.empty() || contig.empty()) {
            cerr << "Error: must provide --sample and --contig OR --seq-id" << endl;
            usage(argv[0]);
            return 1;
        }
        cerr << "Warning: Path name lookup not fully implemented. Using --seq-id directly is recommended." << endl;
        return 1;
    }
    
    if (seq_id >= gbz.index.sequences()) {
        cerr << "Error: sequence ID " << seq_id << " is out of range (max: " << (gbz.index.sequences() - 1) << ")" << endl;
        return 1;
    }
    
    // Load r-index
    FastLocate r_index;
    {
        cerr << "Loading r-index: " << r_index_file << "..." << endl;
        ifstream rin(r_index_file, ios::binary);
        if (!rin) { 
            cerr << "Cannot open r-index: " << r_index_file << endl; 
            return 1; 
        }
        r_index.load_encoded(rin);
    }
    cerr << "R-index loaded successfully. BWT size: " << r_index.bwt_size() << endl;
    
    // Load sampled tag array
    SampledTagArray sampled;
    {
        cerr << "Loading sampled tag array: " << sampled_tags_file << "..." << endl;
        ifstream sin(sampled_tags_file, ios::binary);
        if (!sin) { 
            cerr << "Cannot open sampled.tags: " << sampled_tags_file << endl; 
            return 1; 
        }
        sampled.load(sin);
    }
    cerr << "Sampled tag array loaded successfully. Total runs: " << sampled.total_runs() << endl;
    
    // Extract path sequence from GBWT
    cerr << "Extracting sequence " << seq_id << " from GBWT..." << endl;
    auto path_nodes = gbz.index.extract(seq_id);
    if (path_nodes.empty()) {
        cerr << "Error: sequence " << seq_id << " is empty" << endl;
        return 1;
    }
    cerr << "Extracted " << path_nodes.size() << " nodes from sequence " << seq_id << endl;
    
    // Convert path interval to GBWT sequence interval
    size_t seq_start = interval.first;
    size_t seq_end = interval.second;
    
    // Build full sequence string from path
    cerr << "Building full sequence from path nodes..." << endl;
    string full_seq;
    for (size_t i = 0; i < path_nodes.size(); ++i) {
        gbwt::short_type node_short = path_nodes[i];
        size_t node_id = gbwt::Node::id(node_short);
        bool is_rev = gbwt::Node::is_reverse(node_short);
        
        handle_t handle = gbz.graph.get_handle(node_id, is_rev);
        string node_seq = gbz.graph.get_sequence(handle);
        full_seq += node_seq;
    }
    cerr << "Full sequence length: " << full_seq.size() << " characters" << endl;
    
    if (seq_end >= full_seq.size()) {
        cerr << "Warning: interval end " << seq_end << " exceeds sequence length " << full_seq.size() << ", truncating" << endl;
        seq_end = full_seq.size() - 1;
    }
    
    // Extract query substring from path interval
    string query_seq = full_seq.substr(seq_start, seq_end - seq_start + 1);
    
    if (query_seq.empty()) {
        cerr << "Error: query sequence is empty" << endl;
        return 1;
    }
    cerr << "Query sequence length: " << query_seq.size() << " characters (interval " << seq_start << ".." << seq_end << ")" << endl;
    
    // Step 1: Convert sequence interval to packed text positions in last space
    // last marks packed text positions (not BWT positions)
    cerr << "Converting sequence interval to packed text positions..." << endl;
    size_t text_pos_i = r_index.pack(seq_id, seq_start);  // Start of interval in packed text space
    size_t text_pos_j = r_index.pack(seq_id, seq_end);    // End of interval in packed text space
    
    cerr << "Sequence interval [" << seq_start << ", " << seq_end << "] maps to text positions [" 
         << text_pos_i << ", " << text_pos_j << "] in last space" << endl;
    
    // Step 2: Find successor position on last at or after j
    // Use last_successor(j) to find the first marked position x at or after j
    r_index.ensure_last_rank();
    r_index.ensure_last_select();
    
    auto successor_result = r_index.last_successor(text_pos_j);
    size_t text_pos_x = successor_result.first;   // The marked text position at or after j
    size_t rank_x = successor_result.second;      // The rank of text_pos_x
    
    // Get run ID: run_id = last_to_run[rank_x - 1] (0-indexed)
    size_t run_id = (rank_x > 0 && rank_x <= r_index.last_to_run.size()) 
        ? r_index.last_to_run[rank_x - 1] : 0;
    
    // Find BWT position ISA[x] at the end of this run
    // Use bwt_end_position_of_run to get the BWT end position of the run
    size_t bwt_pos_x = r_index.bwt_end_position_of_run(run_id);
    
    cerr << "Found marked text position x=" << text_pos_x << " (rank=" << rank_x 
         << ", run_id=" << run_id << ", BWT_pos=" << bwt_pos_x << ")" << endl;
    
    // Step 3: Iterate backwards: (x-1, LF(ISA[x])) until x reaches i
    // Collect tags from TAG[ISA[x]] for each x in [i, j]
    cerr << "Iterating backwards through text positions [" << text_pos_i << ", " << text_pos_x << "]..." << endl;
    sampled.ensure_run_rank();
    sampled.ensure_run_select();
    
    unordered_map<uint64_t, TagResult> results;
    size_t current_text_pos = text_pos_x;
    size_t current_bwt_pos = bwt_pos_x;
    
    // Cache sequence lengths
    unordered_map<size_t, size_t> seq_length_cache;
    auto get_seq_length = [&](size_t seq_id) -> size_t {
        auto it = seq_length_cache.find(seq_id);
        if (it != seq_length_cache.end()) {
            return it->second;
        }
        auto path_nodes_seq = gbz.index.extract(seq_id);
        size_t length = 0;
        for (size_t i = 0; i < path_nodes_seq.size(); ++i) {
            gbwt::short_type node_short = path_nodes_seq[i];
            size_t node_id = gbwt::Node::id(node_short);
            bool is_rev = gbwt::Node::is_reverse(node_short);
            handle_t handle = gbz.graph.get_handle(node_id, is_rev);
            length += gbz.graph.get_length(handle);
        }
        seq_length_cache[seq_id] = length;
        return length;
    };
    
    // Iterate backwards from x down to i
    // Only process positions in [i, j] (x might be > j since it's the successor)
    while (current_text_pos >= text_pos_i) {
        // Only process positions in our interval [i, j]
        if (current_text_pos <= text_pos_j) {
            // Get tag for current BWT position
            size_t sampled_run_id = sampled.run_id_at(current_bwt_pos);
            uint64_t tag_val = sampled.run_value(sampled_run_id);
            

            // Do nothing if tag_val == 0
            // If we found a tag, locate the BWT position to get (seq_id, offset)
            if (tag_val != 0) {
                // Locate current BWT position
                size_t rindex_run_id_pos = 0;
                size_t run_start_pos = 0;
                r_index.run_id_and_offset_at(current_bwt_pos, rindex_run_id_pos, run_start_pos);
                size_t packed_pos = r_index.getSample(rindex_run_id_pos);
                
                // Navigate from run start to current_bwt_pos
                for (size_t p = run_start_pos; p < current_bwt_pos; ++p) {
                    packed_pos = r_index.locateNext(packed_pos);
                }
                
                // Unpack to get (seq_id, offset)
                auto pr = r_index.unpack(packed_pos);
                size_t seq_id_found = pr.first;
                size_t offset_from_start = pr.second; // offset is already distance from START
                
                // Unpack current_text_pos to verify it matches
                size_t text_seq_id = r_index.seqId(current_text_pos);
                size_t text_offset = r_index.seqOffset(current_text_pos);
                
                // Calculate offset relative to query start (only for query sequence)
                if (seq_id_found == seq_id && offset_from_start >= seq_start && offset_from_start <= seq_end) {
                    size_t offset_in_query = offset_from_start - seq_start;
                    results[tag_val].query_offsets.push_back(offset_in_query);
                }
                
                cerr << "  Text pos " << current_text_pos << " (seq_id=" << text_seq_id 
                     << ", offset=" << text_offset << ") -> BWT pos " << current_bwt_pos 
                     << " -> (seq_id=" << seq_id_found << ", offset_from_start=" << offset_from_start 
                     << ") -> tag=" << tag_val << endl;
            }
        }
        
        // Move to previous text position: x--
        if (current_text_pos > text_pos_i) {
            current_text_pos--;
            
            // Apply LF-mapping: ISA[x-1] = LF(ISA[x])
            current_bwt_pos = r_index.LF(current_bwt_pos);
        } else {
            break; // Reached i
        }
    }
    
    cerr << "Found " << results.size() << " unique tags from backward iteration" << endl;
    
    if (results.empty()) {
        cout << "seq_id=" << seq_id << "\tinterval=" << seq_start << ".." << seq_end << "\tno_tags" << endl;
        return 0;
    }
    
    // Step 3: For each tag, find all occurrences using wt_gmr select queries
    cerr << "Enumerating all occurrences for " << results.size() << " tags..." << endl;
    const auto& wm = sampled.values();
    
    // Note: seq_length_cache/get_seq_length already defined above
 
    // Cache sequences to avoid recomputation
    unordered_map<size_t, string> seq_cache;
    auto get_sequence = [&](size_t seq_id) -> const string& {
        auto it = seq_cache.find(seq_id);
        if (it != seq_cache.end()) {
            return it->second;
        }
        auto path_nodes_seq = gbz.index.extract(seq_id);
        string full_seq;
        for (size_t i = 0; i < path_nodes_seq.size(); ++i) {
            gbwt::short_type node_short = path_nodes_seq[i];
            size_t node_id = gbwt::Node::id(node_short);
            bool is_rev = gbwt::Node::is_reverse(node_short);
            handle_t handle = gbz.graph.get_handle(node_id, is_rev);
            string node_seq = gbz.graph.get_sequence(handle);
            full_seq += node_seq;
        }
        seq_cache[seq_id] = full_seq;
        return seq_cache[seq_id];
    };
    
    for (auto& kv : results) {
        uint64_t tag_code = kv.first;
        TagResult& res = kv.second;
        
        // Deduplicate query offsets
        sort(res.query_offsets.begin(), res.query_offsets.end());
        res.query_offsets.erase(unique(res.query_offsets.begin(), res.query_offsets.end()), res.query_offsets.end());
        
        // Get total occurrences of this tag using wt_gmr rank
        size_t total_occ = wm.rank(wm.size(), tag_code);
        if (total_occ == 0) {
            cerr << "Warning: Tag " << tag_code << " has 0 occurrences in wavelet matrix" << endl;
            continue;
        }
        
        cerr << "  Tag " << tag_code << ": Found " << total_occ << " tag runs containing this tag" << endl;
        
        // Step 4: For each run containing this tag, locate all (seq_id, offset) pairs
        // Use select queries to find all runs containing this tag
        // We want ALL sequences that have this tag, not just those in the query interval
        for (size_t j = 1; j <= total_occ; ++j) {
            size_t occ_run = wm.select(j, tag_code);  // This is a TAG run_id, not a BWT run_id
            auto run_span = sampled.run_span(occ_run);
            
            cerr << "    Tag Run " << j << "/" << total_occ << ": tag_run_id=" << occ_run 
                 << ", BWT interval=[" << run_span.first << ", " << run_span.second 
                 << "] (size=" << (run_span.second - run_span.first + 1) << ")" << endl;
            
            // Process the entire run to get all sequences with this tag
            // Convert from BWT positions to (seq_id, offset) for all positions in the run
            size_t locate_start = run_span.first;
            size_t locate_end = run_span.second;
            
            // Use manual locating with verification against decompressSA
            // Load decompressSA for verification only (not for actual locating)
            static bool sa_loaded = false;
            static vector<gbwt::size_type> sa;
            if (!sa_loaded) {
                cerr << "      Loading decompressSA for verification..." << endl;
                sa = r_index.decompressSA();
                sa_loaded = true;
            }
            
            // Process all BWT positions in this tag run interval using manual locating
            unordered_set<unsigned long long> seen;
            size_t positions_processed = 0;
            size_t unique_sequences_found = 0;
            size_t mismatches = 0;
            
            // Process positions by r-index runs to optimize
            size_t current_pos = locate_start;
            
            while (current_pos <= locate_end) {
                // Use run_id_and_offset_at to find the run containing current_pos
                // This is the same approach used by locate_encoded
                size_t rindex_run_id = 0;
                size_t run_start_pos = 0;
                r_index.run_id_and_offset_at(current_pos, rindex_run_id, run_start_pos);
                
                // Get sample for this run (packed position at run start)
                size_t packed_pos = r_index.getSample(rindex_run_id);
                
                // Find the end position of this run (start of next run, or end of BWT)
                // We can find this by checking where the next run starts
                size_t rank_at_pos = r_index.last_rank_1(current_pos + 1);
                size_t run_end_pos;
                if (rank_at_pos < r_index.last_to_run.size()) {
                    run_end_pos = r_index.last_select_1(rank_at_pos + 1);
                } else {
                    run_end_pos = r_index.bwt_size() - 1;
                }
                
                // Process positions in this run segment
                size_t segment_start = std::max(current_pos, locate_start);
                size_t segment_end = std::min(run_end_pos, locate_end);
                
                // Navigate from run start to segment_start
                // packed_pos currently points to run_start_pos, we need to advance to segment_start
                for (size_t p = run_start_pos; p < segment_start; ++p) {
                    packed_pos = r_index.locateNext(packed_pos);
                }
                
                // Process all positions in this segment
                for (size_t pos = segment_start; pos <= segment_end; ++pos) {
                    if (pos > segment_start) {
                        packed_pos = r_index.locateNext(packed_pos);
                    }
                    
                    auto pr = r_index.unpack(packed_pos);
                    size_t seq_id_found = pr.first;
                    size_t offset_from_start = pr.second; // offset is already distance from START
                    
                    // Verify against decompressSA
                    bool verified = false;
                    if (pos < sa.size()) {
                        gbwt::size_type sa_packed = sa[pos];
                        auto sa_pr = r_index.unpack(sa_packed);
                        size_t sa_seq_id = sa_pr.first;
                        size_t sa_offset_from_start = sa_pr.second;
                        
                        if (seq_id_found != sa_seq_id || offset_from_start != sa_offset_from_start) {
                            mismatches++;
                            cerr << "        ERROR: BWT[" << pos << "] mismatch! Manual: (seq_id=" << seq_id_found 
                                 << ", offset_from_start=" << offset_from_start 
                                 << "), decompressSA: (seq_id=" << sa_seq_id 
                                 << ", offset_from_start=" << sa_offset_from_start << ")" << endl;
                            // Use decompressSA values since they're correct
                            seq_id_found = sa_seq_id;
                            offset_from_start = sa_offset_from_start;
                        } else {
                            verified = true;
                        }
                    }
                    
                    unsigned long long key = (static_cast<unsigned long long>(seq_id_found) << 32) | static_cast<unsigned long long>(offset_from_start);
                    
                    if (seen.insert(key).second) {
                        res.sequence_offsets[seq_id_found].push_back(offset_from_start);
                        unique_sequences_found++;
                        
                        // Print first few conversions for debugging
                        if (positions_processed < 5) {
                            cerr << "        BWT[" << pos << "] -> (seq_id=" << seq_id_found 
                                 << ", offset_from_start=" << offset_from_start;
                            if (verified) {
                                cerr << ", verified)";
                            } else {
                                cerr << ", corrected from decompressSA)";
                            }
                            cerr << endl;
                        }
                    }
                    positions_processed++;
                }
                
                // Move to next run
                current_pos = segment_end + 1;
            }
            
            if (mismatches > 0) {
                cerr << "      WARNING: Found " << mismatches << " mismatches out of " 
                     << positions_processed << " positions. Used decompressSA values for corrections." << endl;
            }
            
            cerr << "      Processed " << positions_processed << " BWT positions, found " 
                 << unique_sequences_found << " unique (seq_id, offset) pairs" << endl;
        }
    }
    
    // Step 5: Output results
    // For each TAG, show the different sequence IDs that pass through those TAGS
    cerr << "Writing results..." << endl;
    
    cout << "Query Sequence:" << endl;
    cout << "  Sequence ID: " << seq_id << endl;
    cout << "  Interval: " << seq_start << ".." << seq_end << endl;
    cout << "  Sequence: " << query_seq << endl;
    // Compute BWT interval for reporting
    {
        auto __range = r_index.count_encoded(query_seq);
        size_t l = __range.first, r = __range.second;
        cout << "  BWT interval: [" << l << ", " << r << "]" << endl;
    }
    cout << "  Number of tags found: " << results.size() << endl;
    cout << endl;
    
    // Print results for each tag
    for (const auto& kv : results) {
        uint64_t tag_code = kv.first;
        const TagResult& res = kv.second;
        
        // Decode tag code to get node_id and is_rev
        uint64_t decoded = tag_code - 1;
        int64_t node_id = (decoded >> 1) + 1;
        bool is_rev = (decoded & 1) != 0;
        
        cout << "Tag Code: " << tag_code << endl;
        cout << "  Tag Details:" << endl;
        cout << "    Node ID: " << node_id << endl;
        cout << "    Is Reverse: " << (is_rev ? "true" : "false") << endl;
        cout << "    Query offsets (relative to query start): ";
        for (size_t i = 0; i < res.query_offsets.size(); ++i) {
            if (i > 0) cout << ", ";
            cout << res.query_offsets[i];
        }
        cout << endl;
        
        // Show sequence IDs that pass through this tag
        cout << "  Sequence IDs passing through this tag: " << res.sequence_offsets.size() << endl;
        for (const auto& seq_kv : res.sequence_offsets) {
            size_t seq_id_found = seq_kv.first;
            const vector<size_t>& offsets = seq_kv.second;
            
            // Get the full sequence for this seq_id
            const string& full_seq = get_sequence(seq_id_found);
            
            cout << "    Sequence ID: " << seq_id_found << " (occurrences: " << offsets.size() << ")" << endl;
            cout << "      Offsets: ";
            for (size_t i = 0; i < offsets.size() && i < 10; ++i) {
                if (i > 0) cout << ", ";
                cout << offsets[i];
            }
            if (offsets.size() > 10) {
                cout << ", ... (and " << (offsets.size() - 10) << " more)";
            }
            cout << endl;
            
            // Print sequence snippets for each offset
            cout << "      Sequences: ";
            for (size_t i = 0; i < offsets.size() && i < 10; ++i) {
                if (i > 0) cout << ", ";
                size_t offset = offsets[i];
                if (offset < full_seq.size()) {
                    size_t snippet_len = std::min<size_t>(10, full_seq.size() - offset);
                    cout << "\"" << full_seq.substr(offset, snippet_len) << "\"";
                } else {
                    cout << "\"<out_of_range>\"";
                }
            }
            if (offsets.size() > 10) {
                cout << ", ...";
            }
            cout << endl;
        }
        cout << endl;
    }
    
    cerr << "Query completed successfully" << endl;
    
    return 0;
}
