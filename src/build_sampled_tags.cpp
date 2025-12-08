#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sdsl/sd_vector.hpp>

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
    
    // TEST: Build sd_vector with multiset=true and test predecessor queries
    std::cerr << "\n=== TEST: sd_vector with multiset=true ===" << std::endl;
    std::vector<uint64_t> test_positions = {0, 4, 8, 8, 10, 12};
    size_t universe_size = 13; // Make universe size larger than max position
    size_t num_ones = test_positions.size(); // 6 positions total
    
    // Build sd_vector with multiset=true
    sdsl::sd_vector_builder builder(universe_size, num_ones, true); // multiset=true
    for (uint64_t pos : test_positions) {
        builder.set(pos);
    }
    sdsl::sd_vector<> test_sd(builder);
    
    std::cerr << "Built sd_vector with positions: ";
    for (size_t i = 0; i < test_positions.size(); ++i) {
        std::cerr << test_positions[i];
        if (i < test_positions.size() - 1) std::cerr << ", ";
    }
    std::cerr << std::endl;
    
    // Test predecessor queries for positions 0-7
    std::cerr << "\nPredecessor queries:" << std::endl;
    for (size_t pos = 0; pos <= 12; ++pos) {
        auto iter = test_sd.predecessor(pos);
        std::cerr << "  predecessor(" << pos << ") -> ";
        std::cerr << "rank=" << iter->first << ", position=" << iter->second << std::endl;
    }
    std::cerr << "=== END TEST ===" << std::endl;
    
    return 0;
}


