#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "pangenome_index/tag_arrays.hpp"
#include "pangenome_index/sampled_tag_array.hpp"

using namespace std;
using namespace panindexer;

int main(int argc, char** argv) {
    if (argc != 3) {
        cerr << "Usage: ./bin/build_sampled_tags <compact_tags.tags> <sampled.tags>" << endl;
        return 1;
    }
    string in_path = argv[1];
    string out_path = argv[2];

    // Load compact tag arrays (sdsl int_vector-based, with bwt intervals)
    TagArray tags;
    {
        ifstream in(in_path, ios::binary);
        if (!in) { cerr << "Cannot open input: " << in_path << endl; return 1; }
        tags.load_compressed_tags_compact(in);
    }
    std::cerr << "Loaded tags" << std::endl;
    std::cerr << "The size of the tags is: " << tags.bwt_size() << std::endl;

    SampledTagArray sampled;

    // Stream runs from TagArray into SampledTagArray builder
    auto enumerator = [&](const function<void(handlegraph::pos_t,uint64_t)>& sink){
        tags.for_each_run_compact([&](handlegraph::pos_t p, uint64_t len){
            sink(p, len);
        });
    };
    std::cerr << "Building sampled tags" << std::endl;
    sampled.build_from_enumerator(enumerator, tags.bwt_size());
    
    // Print all sampled tags in order (including 0s for gaps)
    cerr << "\nPrinting all sampled tags in order..." << endl;
    sampled.ensure_run_select();
    
    size_t total_runs = sampled.total_runs();
    cerr << "Total runs: " << total_runs << endl;
    
    for (size_t run_id = 0; run_id < total_runs; ++run_id) {
        uint64_t tag_val = sampled.run_value(run_id);
        auto span = sampled.run_span(run_id);
        
        if (tag_val == 0) {
            cerr << "Run " << run_id << ": tag=0 (gap), BWT_pos=[" << span.first << ", " << span.second << "]" << endl;
        } else {
            // Decode tag value: tag_val = 1 + ((node_id - 1) << 1) | is_rev
            uint64_t decoded = tag_val - 1;
            int64_t node_id = (decoded >> 1) + 1;
            bool is_rev = (decoded & 1) != 0;
            
            cerr << "Run " << run_id << ": tag=" << tag_val 
                 << ", node_id=" << node_id 
                 << ", is_rev=" << (is_rev ? "true" : "false")
                 << ", BWT_pos=[" << span.first << ", " << span.second << "]" << endl;
        }
    }
    
    cerr << "Printed " << total_runs << " sampled tags" << endl;
    
    // On macOS, add a small delay to ensure all destructors are called
    // This helps identify if the issue is during destruction
    std::cerr << "About to serialize..." << std::endl;

    // Serialize sampled structure
    ofstream out(out_path, ios::binary);
    if (!out) { cerr << "Cannot open output: " << out_path << endl; return 1; }
    std::cerr << "Calling serialize..." << std::endl;
    sampled.serialize(out);
    std::cerr << "Serialize completed, closing file..." << std::endl;
    out.close();
    std::cerr << "File closed successfully" << std::endl;
    return 0;
}


