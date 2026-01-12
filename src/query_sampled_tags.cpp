#include "pangenome_index/r-index.hpp"
#include "pangenome_index/sampled_tag_array.hpp"
#include <sdsl/wavelet_trees.hpp>
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
using namespace panindexer;
using namespace std::chrono;

static bool debug = true;

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
    cerr << "Usage: " << prog << " <r_index.ri> <sampled.tags> --seq-id ID --interval START..END" << endl;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }
    
    string r_index_file = argv[1];
    string sampled_tags_file = argv[2];
    
    // Parse optional arguments
    size_t seq_id = numeric_limits<size_t>::max();
    pair<size_t, size_t> interval = {0, 0};
    bool has_interval = false;
    
    for (int i = 3; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--interval" && i + 1 < argc) {
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
    
    if (seq_id == numeric_limits<size_t>::max()) {
        cerr << "Error: --seq-id is required" << endl;
        usage(argv[0]);
        return 1;
    }
    
    if (!has_interval) {
        cerr << "Error: --interval is required" << endl;
        usage(argv[0]);
        return 1;
    }
    
    // Timing variables
    auto start_total = high_resolution_clock::now();
    auto start_load_rindex = high_resolution_clock::now();
    
    // Load r-index
    FastLocate r_index;
    {
        if (debug) cerr << "Loading r-index: " << r_index_file << "..." << endl;
        ifstream rin(r_index_file, ios::binary);
        if (!rin) { 
            cerr << "Cannot open r-index: " << r_index_file << endl; 
            return 1; 
        }
        r_index.load_encoded(rin);
    }
    auto end_load_rindex = high_resolution_clock::now();
    auto duration_load_rindex = duration_cast<milliseconds>(end_load_rindex - start_load_rindex);
    if (debug) cerr << "R-index loaded successfully. BWT size: " << r_index.bwt_size() << endl;
    
    auto start_load_tags = high_resolution_clock::now();
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
    auto end_load_tags = high_resolution_clock::now();
    auto duration_load_tags = duration_cast<milliseconds>(end_load_tags - start_load_tags);
    if (debug) cerr << "Sampled tag array loaded successfully. Total runs: " << sampled.total_runs() << endl;
    
    // Convert path interval to sequence interval
    size_t seq_start = interval.first;
    size_t seq_end = interval.second;
    
    // Step 1: Convert sequence interval to packed text positions in last space
    // last marks packed text positions (not BWT positions)
    if (debug) cerr << "Converting sequence interval to packed text positions..." << endl;
    size_t text_pos_i = r_index.pack(seq_id, seq_start);  // Start of interval in packed text space
    size_t text_pos_j = r_index.pack(seq_id, seq_end);    // End of interval in packed text space
    
    if (debug) cerr << "Sequence interval [" << seq_start << ", " << seq_end << "] maps to text positions [" 
         << text_pos_i << ", " << text_pos_j << "] in last space" << endl;
    
    // Step 2: Find successor position on last at or after j
    // Use last_successor(j) to find the first marked position x at or after j
    r_index.ensure_last_rank();
    r_index.ensure_last_select();
    
    auto successor_result = r_index.last_successor(text_pos_j);
    size_t text_pos_x = successor_result.first;   // The marked text position at or after j
    size_t rank_x = successor_result.second;      // 0-based index into last_to_run
    
    // Get run ID: rank_x is already 0-based index into last_to_run
    size_t run_id = (rank_x < r_index.last_to_run.size()) 
        ? r_index.last_to_run[rank_x] : 0;
    
    // Find BWT position ISA[x] at the end of this run
    // Use bwt_end_position_of_run to get the BWT end position of the run
    size_t bwt_pos_x = r_index.bwt_end_position_of_run(run_id);
    
    if (debug) cerr << "Found marked text position x=" << text_pos_x << " (rank=" << rank_x 
         << ", run_id=" << run_id << ", BWT_pos=" << bwt_pos_x << ")" << endl;
    
    // Check if text_pos_x is after the sequence end, and if so, start from BWT[seq_id]
    size_t bwt_seq_end = seq_id;  // BWT position of $ at end of sequence seq_id
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
    
    if (debug) cerr << "Sequence end check: text_pos_x=" << text_pos_x 
         << ", text_pos_seq_end=" << text_pos_seq_end << " (BWT[" << bwt_seq_end << "])" << endl;
    
    // Check if text_pos_x (BWT rank from last array) is greater than sequence end
    if (text_pos_x > text_pos_seq_end) {
        cerr << "Error: BWT rank calculated from last array (text_pos_x=" << text_pos_x 
             << ") is greater than sequence end (text_pos_seq_end=" << text_pos_seq_end << ")" << endl;
        return 1;
    }
    
    // Step 3: Iterate backwards: (x-1, LF(ISA[x])) until x reaches i
    // Collect tags from TAG[ISA[x]] for each x in [i, j]
    auto start_find_tags = high_resolution_clock::now();
    if (debug) cerr << "Iterating backwards through text positions [" << text_pos_i << ", " << text_pos_x << "]..." << endl;
    sampled.ensure_run_rank();
    sampled.ensure_run_select();
    
    unordered_map<uint64_t, TagResult> results;
    size_t current_text_pos = text_pos_x;
    size_t current_bwt_pos = bwt_pos_x;
    size_t positions_checked = 0;
    size_t positions_with_tags = 0;
    
    // Iterate backwards from x down to i
    // Only process positions in [i, j] (x might be > j since it's the successor)
    while (current_text_pos >= text_pos_i) {
        // Only process positions in our interval [i, j]
        if (current_text_pos <= text_pos_j) {
            positions_checked++;
            // Get tag for current BWT position
            size_t sampled_run_id = sampled.run_id_at(current_bwt_pos);
            uint64_t tag_val = sampled.run_value(sampled_run_id);
            

            // Do nothing if tag_val == 0
            // If we found a tag, unpack current_text_pos to get (seq_id, offset)
            if (tag_val != 0) {
                positions_with_tags++;
                // Unpack current_text_pos directly (it's already a packed position)
                size_t seq_id_found = r_index.seqId(current_text_pos);
                size_t offset_from_start = r_index.seqOffset(current_text_pos);
                
                // Calculate offset relative to query start (only for query sequence)
                if (seq_id_found == seq_id && offset_from_start >= seq_start && offset_from_start <= seq_end) {
                    size_t offset_in_query = offset_from_start - seq_start;
                    results[tag_val].query_offsets.push_back(offset_in_query);
                }
                
                if (debug) cerr << "  Text pos " << current_text_pos << " (seq_id=" << seq_id_found 
                     << ", offset=" << offset_from_start << ") -> BWT pos " << current_bwt_pos 
                     << " -> tag=" << tag_val << endl;
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
    auto end_find_tags = high_resolution_clock::now();
    auto duration_find_tags = duration_cast<milliseconds>(end_find_tags - start_find_tags);
    if (debug) cerr << "Found " << results.size() << " unique tags from backward iteration" << endl;
    
    if (results.empty()) {
        auto end_total = high_resolution_clock::now();
        auto duration_total = duration_cast<milliseconds>(end_total - start_total);
        cout << "seq_id=" << seq_id << "\tinterval=" << seq_start << ".." << seq_end << "\tno_tags" << endl;
        cerr << "\n=== Query Statistics ===" << endl;
        cerr << "  Sequence ID: " << seq_id << endl;
        cerr << "  Interval: " << seq_start << ".." << seq_end << " (length: " << (seq_end - seq_start + 1) << ")" << endl;
        cerr << "  Positions checked: " << positions_checked << endl;
        cerr << "  Positions with tags: " << positions_with_tags << endl;
        cerr << "  Unique tags found: 0" << endl;
        cerr << "\n=== Timing Statistics ===" << endl;
        cerr << "  Load r-index: " << duration_load_rindex.count() << " ms" << endl;
        cerr << "  Load sampled tags: " << duration_load_tags.count() << " ms" << endl;
        cerr << "  Find tags in interval: " << duration_find_tags.count() << " ms" << endl;
        cerr << "  Total time: " << duration_total.count() << " ms" << endl;
        return 0;
    }
    
    // Step 3: For each tag, find all occurrences using wt_gmr select queries
    auto start_find_sequences = high_resolution_clock::now();
    if (debug) cerr << "Enumerating all occurrences for " << results.size() << " tags..." << endl;
    const auto& wm = sampled.values();

    
    size_t total_tag_runs = 0;
    size_t total_bwt_positions_processed = 0;
    size_t total_unique_sequences = 0;
    size_t total_occurrences = 0;
    
    for (auto& kv : results) {
        uint64_t tag_code = kv.first;
        TagResult& res = kv.second;
        
        // Deduplicate query offsets
        sort(res.query_offsets.begin(), res.query_offsets.end());
        res.query_offsets.erase(unique(res.query_offsets.begin(), res.query_offsets.end()), res.query_offsets.end());
        
        // Get total occurrences of this tag using wt_gmr rank
        size_t total_occ = wm.rank(wm.size(), tag_code);
        if (total_occ == 0) {
            if (debug) cerr << "Warning: Tag " << tag_code << " has 0 occurrences in wavelet matrix" << endl;
            continue;
        }
        
        total_tag_runs += total_occ;
        if (debug) cerr << "  Tag " << tag_code << ": Found " << total_occ << " tag runs containing this tag" << endl;
        
        // Step 4: For each run containing this tag, locate all (seq_id, offset) pairs
        // Use select queries to find all runs containing this tag
        // We want ALL sequences that have this tag, not just those in the query interval
        bool first_run_is_gap = sampled.is_first_run_gap();
        for (size_t j = 1; j <= total_occ; ++j) {
            size_t tag_run_id = wm.select(j, tag_code);  // This is a TAG run_id (0-indexed into wt_gmr, non-gap only)
            // Convert TAG run_id to BWT run_id (which includes gaps)
            // If first_run_is_gap=false: BWT run_id = 2 * TAG run_id (even = normal tags)
            // If first_run_is_gap=true: BWT run_id = 2 * TAG run_id + 1 (odd = normal tags)
            size_t bwt_run_id = 2 * tag_run_id + 1 - static_cast<size_t>(first_run_is_gap);
            auto run_span = sampled.run_span(bwt_run_id);
            
            if (debug) cerr << "    Tag Run " << j << "/" << total_occ << ": tag_run_id=" << tag_run_id 
                 << ", bwt_run_id=" << bwt_run_id
                 << ", BWT interval=[" << run_span.first << ", " << run_span.second 
                 << "] (size=" << (run_span.second - run_span.first + 1) << ")" << endl;
            
            // Process the entire run to get all sequences with this tag
            size_t locate_start = run_span.first;
            size_t locate_end = run_span.second;
            
            size_t positions_processed = 0;
            size_t unique_sequences_found = 0;

            size_t packed_pos;
            {
                size_t rindex_run_id = 0;
                size_t run_start_pos = 0;
                r_index.run_id_and_offset_at(locate_start, rindex_run_id, run_start_pos);
                packed_pos = r_index.getSample(rindex_run_id);
                
                // Navigate from run start to locate_start
                for (size_t p = run_start_pos; p < locate_start; ++p) {
                    packed_pos = r_index.locateNext(packed_pos);
                }
            }
            
            // Process all positions sequentially
            for (size_t pos = locate_start; pos <= locate_end; ++pos) {
                if (pos > locate_start) {
                    packed_pos = r_index.locateNext(packed_pos);
                }
                
                // Unpack to get (seq_id, offset)
                auto pr = r_index.unpack(packed_pos);
                size_t seq_id_found = pr.first;
                size_t offset_from_start = pr.second; // offset is already distance from START
                
                // Each BWT position maps to a unique (seq_id, offset) pair
                res.sequence_offsets[seq_id_found].push_back(offset_from_start);
                unique_sequences_found++;
                
                // Print first few conversions for debugging
                if (debug && positions_processed < 5) {
                    cerr << "        BWT[" << pos << "] -> (seq_id=" << seq_id_found 
                         << ", offset_from_start=" << offset_from_start << ")" << endl;
                }
                positions_processed++;
            }
            
            total_bwt_positions_processed += positions_processed;
            total_unique_sequences += unique_sequences_found;
            total_occurrences += positions_processed;
            
            if (debug) cerr << "      Processed " << positions_processed << " BWT positions, found " 
                 << unique_sequences_found << " unique (seq_id, offset) pairs" << endl;
        }
    }
    auto end_find_sequences = high_resolution_clock::now();
    auto duration_find_sequences = duration_cast<milliseconds>(end_find_sequences - start_find_sequences);
    
    // Step 5: Output results
    // For each TAG, show the different sequence IDs that pass through those TAGS
    if (debug) cerr << "Writing results..." << endl;
    
    cout << "Query Sequence:" << endl;
    cout << "  Sequence ID: " << seq_id << endl;
    cout << "  Interval: " << seq_start << ".." << seq_end << endl;
    cout << "  Number of tags found: " << results.size() << endl;
    cout << endl;
    
    // Sort tags by tag_code before printing
    vector<pair<uint64_t, TagResult>> sorted_results(results.begin(), results.end());
    sort(sorted_results.begin(), sorted_results.end(), 
         [](const pair<uint64_t, TagResult>& a, const pair<uint64_t, TagResult>& b) {
             return a.first < b.first;
         });
    
    // Print results for each tag
    for (const auto& kv : sorted_results) {
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
        }
        cout << endl;
    }
    
    auto end_total = high_resolution_clock::now();
    auto duration_total = duration_cast<milliseconds>(end_total - start_total);
    
    // Print statistics
    cerr << "\n=== Query Statistics ===" << endl;
    cerr << "  Sequence ID: " << seq_id << endl;
    cerr << "  Interval: " << seq_start << ".." << seq_end << " (length: " << (seq_end - seq_start + 1) << ")" << endl;
    cerr << "  Positions checked: " << positions_checked << endl;
    cerr << "  Positions with tags: " << positions_with_tags << endl;
    cerr << "  Unique tags found: " << results.size() << endl;
    cerr << "  Total tag runs: " << total_tag_runs << endl;
    cerr << "  Total BWT positions processed: " << total_bwt_positions_processed << endl;
    cerr << "  Total unique sequences: " << total_unique_sequences << endl;
    cerr << "  Total occurrences: " << total_occurrences << endl;
    if (results.size() > 0) {
        cerr << "  Average occurrences per tag: " << fixed << setprecision(2) 
             << (double)total_occurrences / results.size() << endl;
        cerr << "  Average sequences per tag: " << fixed << setprecision(2) 
             << (double)total_unique_sequences / results.size() << endl;
    }
    
    cerr << "\n=== Timing Statistics ===" << endl;
    cerr << "  Load r-index: " << duration_load_rindex.count() << " ms" << endl;
    cerr << "  Load sampled tags: " << duration_load_tags.count() << " ms" << endl;
    cerr << "  Find tags in interval: " << duration_find_tags.count() << " ms" << endl;
    cerr << "  Find all sequences for tags: " << duration_find_sequences.count() << " ms" << endl;
    cerr << "  Total time: " << duration_total.count() << " ms" << endl;
    
    if (debug) cerr << "Query completed successfully" << endl;
    
    return 0;
}
