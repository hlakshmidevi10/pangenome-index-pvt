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

using namespace std;
using namespace panindexer;
using namespace gbwtgraph;

struct TagResult {
    vector<size_t> query_offsets; // offsets within the query interval (relative to l)
    vector<pair<size_t,size_t>> aligned_positions; // (seq_id, seq_offset)
    vector<pair<size_t,size_t>> other_aligned_positions; // occurrences from runs not overlapping [l, r]
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

// Find GBWT sequence id from path name components
static size_t find_path_sequence_id(const GBZ& gbz, const string& sample, const string& contig, size_t haplotype) {
    // Build path name: sample#contig#haplotype or similar format
    // Note: GBWT path naming convention may vary; adjust based on actual format
    string path_name = sample + "#" + contig;
    if (haplotype > 0) {
        path_name += "#" + to_string(haplotype);
    }
    
    // Search through GBWT paths
    for (size_t i = 0; i < gbz.index.sequences(); ++i) {
        // Try to match path name (exact matching depends on GBWT metadata)
        // For now, we'll search by extracting paths and checking
        auto path_nodes = gbz.index.extract(i);
        if (path_nodes.empty()) continue;
        
        // Check if this matches (simplified - may need adjustment based on actual path naming)
        // This is a placeholder - actual implementation depends on GBWT metadata structure
    }
    
    // Fallback: if not found, return 0 (will need error handling)
    return 0;
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
        // Find sequence ID from path name (simplified - may need adjustment)
        // For now, we'll iterate through sequences to find a match
        // TODO: Implement proper path name lookup if GBWT metadata supports it
        cerr << "Warning: Path name lookup not fully implemented. Using --seq-id directly is recommended." << endl;
        // For demonstration, we'll assume seq_id needs to be provided
        return 1;
    }
    
    if (seq_id >= gbz.index.sequences()) {
        cerr << "Error: sequence ID " << seq_id << " is out of range (max: " << (gbz.index.sequences() - 1) << ")" << endl;
        return 1;
    }
    
    // Load r-index (encoded or legacy supported by same API)
    FastLocate r_index;
    {
        cerr << "Loading r-index: " << r_index_file << "..." << endl;
        ifstream rin(r_index_file, ios::binary);
        if (!rin) { cerr << "Cannot open r-index: " << r_index_file << endl; return 1; }
        r_index.load_encoded(rin);
    }
    cerr << "R-index loaded successfully. BWT size: " << r_index.bwt_size() << endl;
    
    // Load sampled tag array
    SampledTagArray sampled;
    {
        cerr << "Loading sampled tag array: " << sampled_tags_file << "..." << endl;
        ifstream sin(sampled_tags_file, ios::binary);
        if (!sin) { cerr << "Cannot open sampled.tags: " << sampled_tags_file << endl; return 1; }
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
    // The interval is relative to sequence positions (characters) in the path
    size_t seq_start = interval.first;
    size_t seq_end = interval.second;
    
    // Build full sequence string from path
    cerr << "Building full sequence from path nodes..." << endl;
    string full_seq;
    vector<size_t> node_start_positions; // Track where each node starts in the sequence
    node_start_positions.push_back(0);
    
    for (size_t i = 0; i < path_nodes.size(); ++i) {
        gbwt::short_type node_short = path_nodes[i];
        size_t node_id = gbwt::Node::id(node_short);
        bool is_rev = gbwt::Node::is_reverse(node_short);
        
        handle_t handle = gbz.graph.get_handle(node_id, is_rev);
        string node_seq = gbz.graph.get_sequence(handle);
        full_seq += node_seq;
        if (i + 1 < path_nodes.size()) {
            node_start_positions.push_back(full_seq.size());
        }
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
    
    // Map query sequence to BWT interval using r-index count
    cerr << "Mapping query sequence to BWT interval..." << endl;
    gbwt::range_type range = r_index.count_encoded(query_seq);
    if (range.first > range.second) {
        cout << "seq_id=" << seq_id << "\tinterval=" << seq_start << ".." << seq_end << "\tno_matches" << endl;
        return 0;
    }
    size_t l = range.first, r = range.second;
    cerr << "BWT interval: [" << l << ", " << r << "] (size: " << (r - l + 1) << ")" << endl;
    
    // Step: successor query on last at or after r, and corresponding BWT position
    cerr << "Finding successor position on last..." << endl;
    r_index.ensure_last_rank();
    r_index.ensure_last_select();
    size_t ones_upto_r = r_index.last_rank_1(r + 1);
    bool r_is_end = false;
    if (r < r_index.last.size()) { 
        r_is_end = (r_index.last[r] == 1); 
    }
    size_t sel_k = r_is_end ? ones_upto_r : (ones_upto_r + 1);
    size_t succ_bwt_pos = (sel_k == 0 ? r : r_index.last_select_1(sel_k));
    cerr << "Successor position: " << succ_bwt_pos << endl;
    
    // Get run_id from last_to_run for the successor position
    // The last_to_run array maps rank of last to run_id
    size_t run_id_start = (sel_k > 0 && sel_k <= r_index.last_to_run.size()) ? 
                          r_index.last_to_run[sel_k - 1] : r_index.last_to_run[ones_upto_r];
    
    // Now iterate backwards with LF-mapping from succ_bwt_pos down to l
    // We'll collect all non-gap tags encountered in the query interval [l, r]
    cerr << "Collecting tags in BWT interval [" << l << ", " << r << "]..." << endl;
    unordered_map<uint64_t, TagResult> results;
    unordered_map<uint64_t, unordered_set<size_t>> present_runs;
    
    // Start from successor position and iterate backwards using LF-mapping
    // LF-mapping: LF(i) = C[c] + rank(i, c) where c = BWT[i]
    // To go backwards, we use psi: psi(i) gives (BWT[i], LF(i))
    // But we need to iterate backwards through the interval, so we'll traverse positions
    
    // Approach: iterate from r backwards to l using BWT characters
    // For each position in [l, r], check the tag at that position
    
    // More efficient: iterate through runs that overlap [l, r]
    size_t rid_l = sampled.run_id_at(l);
    size_t rid_r = sampled.run_id_at(r);
    cerr << "Run IDs covering interval: [" << rid_l << ", " << rid_r << "]" << endl;
    
    // First pass: collect tags from runs overlapping [l, r]
    size_t tags_found = 0;
    for (size_t run_id = rid_l; run_id <= rid_r; ++run_id) {
        uint64_t val = sampled.run_value(run_id);
        if (val == 0) continue; // skip gaps
        auto span = sampled.run_span(run_id);
        size_t s = span.first, e = span.second;
        size_t a = max(s, l);
        size_t b = min(e, r);
        if (a <= b) {
            size_t offset_in_query = a - l;
            results[val].query_offsets.push_back(offset_in_query);
            present_runs[val].insert(run_id);
            tags_found++;
        }
    }
    cerr << "Found " << tags_found << " tag occurrences in " << results.size() << " unique tags" << endl;
    cerr.flush();
    
    // NOTE: The forward pass above already collected tags from runs overlapping [l, r]
    // The backward iteration below is supposed to use LF-mapping, but for now we'll
    // skip it if tags were already found, to isolate the hanging issue
    // TODO: Implement proper backward LF-mapping iteration
    
    cerr << "About to start backward iteration..." << endl;
    cerr.flush();
    
    // Explicitly initialize rank/select supports before backward iteration
    cerr << "Initializing sampled tag array supports..." << endl;
    cerr.flush();
    sampled.ensure_run_rank();
    sampled.ensure_run_select();
    cerr << "Supports initialized successfully" << endl;
    cerr.flush();
    
    // Iterate backwards from succ_bwt_pos down to l using LF-mapping
    // We start at the successor position (first marked position at or after r)
    // and work backwards through the interval [l, r], collecting non-gap tags
    // This follows the backward LF-mapping traversal as specified
    
    // Start from r (or successor if it's within our interval)
    // The successor position might be slightly after r, so we start from r
    cerr << "Computing current_pos..." << endl;
    cerr.flush();
    size_t current_pos = min(succ_bwt_pos, r);
    cerr << "current_pos = " << current_pos << endl;
    cerr.flush();
    if (current_pos > r) {
        // If successor is after r, we need to go back to r
        current_pos = r;
    }
    
    // Track which runs we've already processed to avoid duplicates
    unordered_set<size_t> processed_runs;
    
    // Debug: Print initial state
    cerr << "Starting backward iteration from position " << current_pos << " down to " << l << "..." << endl;
    cerr.flush();
    
    {
        cerr << "  Getting initial run_id..." << endl;
        cerr.flush();
        size_t initial_run_id = sampled.run_id_at(current_pos);
        cerr << "  Got initial run_id: " << initial_run_id << endl;
        cerr.flush();
        
        cerr << "  Getting initial run_span..." << endl;
        cerr.flush();
        auto initial_span = sampled.run_span(initial_run_id);
        cerr << "  Got initial run_span: [" << initial_span.first << ", " << initial_span.second << "]" << endl;
        cerr.flush();
        
        cerr << "  Getting initial tag value..." << endl;
        cerr.flush();
        uint64_t initial_tag_val = sampled.run_value(initial_run_id);
        cerr << "  Initial state: run_id=" << initial_run_id 
             << ", run_span=[" << initial_span.first << ", " << initial_span.second << "]"
             << ", tag_val=" << initial_tag_val << endl;
        cerr.flush();
    }
    
    size_t iteration_count = 0;
    const size_t MAX_ITERATIONS = 10000; // Safety limit
    
    while (current_pos >= l && current_pos <= r_index.bwt_size() && iteration_count < MAX_ITERATIONS) {
        iteration_count++;
        
        // Print every iteration for debugging when there are few tags
        if (iteration_count <= 10 || iteration_count % 100 == 0) {
            cerr << "  Iteration " << iteration_count << ": position=" << current_pos << ", run_id=" << sampled.run_id_at(current_pos) << endl;
        }
        
        // Get run_id for current position
        size_t run_id = sampled.run_id_at(current_pos);
        
        // Skip if we've already processed this run
        if (processed_runs.find(run_id) != processed_runs.end()) {
            cerr << "    Already processed run " << run_id << ", moving to previous run" << endl;
            // Move to start of previous run
            if (run_id > 0) {
                auto prev_span = sampled.run_span(run_id - 1);
                current_pos = prev_span.second;
                cerr << "    Moved to previous run end position: " << current_pos << endl;
                if (current_pos < l) {
                    cerr << "  Reached start of interval (previous run ends at " << current_pos << " < " << l << ")" << endl;
                    break;
                }
            } else {
                cerr << "  Reached first run (run_id = 0)" << endl;
                break;
            }
            continue;
        }
        processed_runs.insert(run_id);
        
        // Get tag value for this run
        uint64_t tag_val = sampled.run_value(run_id);
        
        // Record non-gap tags that overlap [l, r]
        if (tag_val != 0) {
            auto span = sampled.run_span(run_id);
            size_t s = span.first, e = span.second;
            size_t overlap_start = max(s, l);
            size_t overlap_end = min(e, r);
            if (overlap_start <= overlap_end) {
                size_t offset_in_query = overlap_start - l;
                results[tag_val].query_offsets.push_back(offset_in_query);
                present_runs[tag_val].insert(run_id);
                cerr << "    Found tag " << tag_val << " at run " << run_id << " (span [" << s << ", " << e << "])" << endl;
            }
        }
        
        // Move backwards: go to start of current run, then to previous run
        auto span = sampled.run_span(run_id);
        if (span.first <= l) {
            // We've reached or passed the start of the interval
            cerr << "  Reached start of interval (run " << run_id << " starts at " << span.first << " <= " << l << ")" << endl;
            break;
        }
        
        // Move to the end of the previous run
        if (run_id > 0) {
            auto prev_span = sampled.run_span(run_id - 1);
            size_t prev_pos = prev_span.second;
            
            cerr << "    Run " << run_id << " span: [" << span.first << ", " << span.second << "]" << endl;
            cerr << "    Previous run " << (run_id - 1) << " span: [" << prev_span.first << ", " << prev_span.second << "]" << endl;
            
            // Safety check: ensure we're actually moving backwards
            if (prev_pos >= current_pos) {
                cerr << "  Warning: Not moving backwards! prev_pos=" << prev_pos << " >= current_pos=" << current_pos << endl;
                cerr << "    Current run_id=" << run_id << ", prev_run_id=" << (run_id - 1) << endl;
                break;
            }
            
            current_pos = prev_pos;
            cerr << "    Moved to position " << current_pos << " (end of run " << (run_id - 1) << ")" << endl;
            
            if (current_pos < l) {
                cerr << "  Previous run ends before interval start (" << current_pos << " < " << l << ")" << endl;
                break;
            }
            // Additional safety check
            if (current_pos >= r_index.bwt_size()) {
                cerr << "  Error: current_pos out of bounds: " << current_pos << " >= " << r_index.bwt_size() << endl;
                break;
            }
        } else {
            cerr << "  Reached first run (run_id = 0)" << endl;
            break;
        }
    }
    
    if (iteration_count >= MAX_ITERATIONS) {
        cerr << "Warning: Backward iteration reached maximum iterations (" << MAX_ITERATIONS << ")" << endl;
    } else {
        cerr << "Backward iteration completed after " << iteration_count << " iterations" << endl;
    }
    
    // For each tag code present, enumerate all occurrences across the full sampled array via wavelet matrix
    if (results.empty()) {
        cerr << "No tags found in interval, skipping enumeration" << endl;
    } else {
        cerr << "Enumerating all occurrences for " << results.size() << " tags..." << endl;
        size_t tag_idx = 0;
        for (auto& kv : results) {
            tag_idx++;
            uint64_t tag_code = kv.first;
            TagResult& res = kv.second;
            
            if (tag_idx % 10 == 0 || tag_idx == 1) {
                cerr << "  Processing tag " << tag_idx << "/" << results.size() << " (code: " << tag_code << ")..." << endl;
            }
            
            // dedup offsets
            sort(res.query_offsets.begin(), res.query_offsets.end());
            res.query_offsets.erase(unique(res.query_offsets.begin(), res.query_offsets.end()), res.query_offsets.end());
            
            const auto& wm = sampled.values();
            // total occurrences of tag_code
            size_t total_occ = wm.rank(wm.size(), tag_code);
            if (total_occ == 0) {
                cerr << "    Tag " << tag_code << " has 0 occurrences in wavelet matrix" << endl;
                continue;
            }
            cerr << "    Tag " << tag_code << " has " << total_occ << " total occurrences" << endl;
            
            // Gather aligned positions for each run where this tag appears
            // We need to locate each position individually to get full (seq_id, offset) pairs
            unordered_set<unsigned long long> seen_all;
            unordered_set<unsigned long long> seen_other;
            for (size_t j = 1; j <= total_occ; ++j) {
                size_t occ_run = wm.select(j, tag_code);
                auto sp = sampled.run_span(occ_run);
                
                // Check if this run overlaps the query interval [l, r]
                auto it = present_runs.find(tag_code);
                bool is_present_run = (it != present_runs.end() && it->second.find(occ_run) != it->second.end());
                
                // Determine which interval to locate: only positions in [l, r] for runs overlapping query
                size_t locate_start, locate_end;
                if (is_present_run) {
                    // For runs overlapping query interval, only locate positions within [l, r]
                    locate_start = max(sp.first, l);
                    locate_end = min(sp.second, r);
                } else {
                    // For runs not overlapping query interval, locate entire run
                    locate_start = sp.first;
                    locate_end = sp.second;
                }
                
                if (locate_start <= locate_end) {
                    // Locate each position individually to get full packed positions (seq_id, offset)
                    // Optimize by caching run info when positions are in the same run
                    size_t current_packed = 0;
                    size_t cached_run_start = 0;
                    size_t cached_pos = locate_start;
                    size_t cached_run_id = numeric_limits<size_t>::max();
                    
                    for (size_t pos = locate_start; pos <= locate_end; pos++) {
                        // Check if we need to recompute (new run or first position)
                        size_t rank = r_index.last_rank_1(pos + 1);
                        size_t run_id = (rank > 0 && rank <= r_index.last_to_run.size()) 
                            ? r_index.last_to_run[rank - 1] : 0;
                        
                        if (run_id != cached_run_id || pos == locate_start) {
                            // New run or first position: recompute from run start
                            cached_run_id = run_id;
                            current_packed = r_index.getSample(run_id);
                            cached_run_start = (rank > 0) ? r_index.last_select_1(rank - 1) + 1 : 0;
                            
                            // Navigate from run_start to pos
                            for (size_t p = cached_run_start; p < pos; p++) {
                                current_packed = r_index.locateNext(current_packed);
                            }
                            cached_pos = pos;
                        } else {
                            // Same run: continue from previous position
                            current_packed = r_index.locateNext(current_packed);
                            cached_pos = pos;
                        }
                        
                        // Now current_packed is the packed position for BWT position pos
                        auto pr = r_index.unpack(current_packed);
                        unsigned long long key = (static_cast<unsigned long long>(pr.first) << 32) | static_cast<unsigned long long>(pr.second);
                        
                        if (seen_all.insert(key).second) {
                            res.aligned_positions.emplace_back(pr.first, pr.second);
                        }
                        // If this occurrence comes from a run not overlapping the current [l, r], treat as "other"
                        if (!is_present_run) {
                            if (seen_other.insert(key).second) {
                                res.other_aligned_positions.emplace_back(pr.first, pr.second);
                            }
                        }
                    }
                }
            }
            cerr << "    Tag " << tag_code << ": " << res.aligned_positions.size() << " aligned positions, " 
                 << res.other_aligned_positions.size() << " other positions" << endl;
        }
        cerr << "Completed enumeration for all tags" << endl;
    }
    
    // Print summary for this path interval
    cerr << "Writing results..." << endl;
    
    // Print the query sequence
    cout << "Query Sequence:" << endl;
    cout << "  Sequence ID: " << seq_id << endl;
    cout << "  Interval: " << seq_start << ".." << seq_end << endl;
    cout << "  Sequence: " << query_seq << endl;
    cout << "  BWT interval: [" << l << ", " << r << "]" << endl;
    cout << "  Number of tags found: " << results.size() << endl;
    cout << endl;
    
    // Print results for each tag
    const auto& wm = sampled.values();
    for (const auto& kv : results) {
        uint64_t code = kv.first;
        const TagResult& res = kv.second;
        
        // Decode tag code to get node_id and is_rev
        // Tag encoding: 1 + ((node_id - 1) << 1) | is_rev
        uint64_t decoded = code - 1;
        int64_t node_id = (decoded >> 1) + 1;
        bool is_rev = (decoded & 1) != 0;
        
        // Get total occurrences of this tag
        size_t total_occ = wm.rank(wm.size(), code);
        
        cout << "Tag Code: " << code << endl;
        cout << "  Tag Details:" << endl;
        cout << "    Node ID: " << node_id << endl;
        cout << "    Is Reverse: " << (is_rev ? "true" : "false") << endl;
        cout << "    Total occurrences in BWT: " << total_occ << endl;
        cout << "    Occurrences in query interval: " << res.query_offsets.size() << endl;
        
        // Show runs that contain this tag and overlap the query interval
        auto it_runs = present_runs.find(code);
        if (it_runs != present_runs.end() && !it_runs->second.empty()) {
            cout << "  Runs overlapping query interval [" << l << ", " << r << "]:" << endl;
            cout << "    Count: " << it_runs->second.size() << endl;
            vector<size_t> run_ids_vec(it_runs->second.begin(), it_runs->second.end());
            sort(run_ids_vec.begin(), run_ids_vec.end());
            for (size_t i = 0; i < run_ids_vec.size(); ++i) {
                size_t run_id = run_ids_vec[i];
                auto span = sampled.run_span(run_id);
                size_t overlap_start = max(span.first, l);
                size_t overlap_end = min(span.second, r);
                cout << "    [" << i << "] Run ID: " << run_id 
                     << ", BWT span: [" << span.first << ", " << span.second << "]"
                     << ", Overlap with query: [" << overlap_start << ", " << overlap_end << "]" << endl;
            }
        }
        
        cout << "  Query offsets (relative to query start): ";
        for (size_t i = 0; i < res.query_offsets.size(); ++i) {
            if (i) cout << ", ";
            cout << res.query_offsets[i];
        }
        cout << endl;
        
        // Print aligned positions (from runs overlapping query interval)
        cout << "  Aligned positions (from runs overlapping query interval):" << endl;
        cout << "    Count: " << res.aligned_positions.size() << endl;
        if (!res.aligned_positions.empty()) {
            for (size_t i = 0; i < res.aligned_positions.size(); ++i) {
                cout << "    [" << i << "] Sequence ID: " << res.aligned_positions[i].first 
                     << ", Offset: " << res.aligned_positions[i].second << endl;
            }
        }
        
        // Print other aligned positions (from runs not overlapping query interval)
        cout << "  Other aligned positions (from runs not overlapping query interval):" << endl;
        cout << "    Count: " << res.other_aligned_positions.size() << endl;
        if (!res.other_aligned_positions.empty()) {
            for (size_t i = 0; i < res.other_aligned_positions.size(); ++i) {
                cout << "    [" << i << "] Sequence ID: " << res.other_aligned_positions[i].first 
                     << ", Offset: " << res.other_aligned_positions[i].second << endl;
            }
        }
        cout << endl;
    }
    cerr << "Query completed successfully" << endl;
    
    return 0;
}
