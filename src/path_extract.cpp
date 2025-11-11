//
// Created by Harieasswar Lakshmidevi on 11/3/25.
//

//
// Created by Harieasswar Lakshmidevi on 11/3/25.
//

#include "pangenome_index/algorithm.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <queue>
#include "pangenome_index/r-index.hpp"

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <gbz_file> <output_file>" << std::endl;
        std::cerr << "Extract all path names from GBZ file and write to output file" << std::endl;
        return 1;
    }

    std::string gbz_file = std::string(argv[1]);
    std::string output_file = std::string(argv[2]);

    // Load the GBZ file
    gbwtgraph::GBZ gbz;
    std::cerr << "Loading GBZ file: " << gbz_file << std::endl;

    try {
        sdsl::simple_sds::load_from(gbz, gbz_file);
    } catch (const std::exception& e) {
        std::cerr << "Error loading GBZ file: " << e.what() << std::endl;
        return 1;
    }

    // Open output file
    std::ofstream output(output_file);
    if (!output.is_open()) {
        std::cerr << "Error: Could not open output file: " << output_file << std::endl;
        return 1;
    }

    std::cerr << "GBZ total sequences: " << gbz.index.sequences() << std::endl;
    std::cerr << "GBZ graph path count: " << gbz.graph.get_path_count() << std::endl;
    std::cerr << "GBZ graph node count: " << gbz.graph.get_node_count() << std::endl;
    std::cerr << "GBZ graph edge count: " << gbz.graph.get_edge_count() << std::endl;

    std::cerr << "Extracting path names..." << std::endl;

    // TODO: check the order and the output of the paths with vg gbwt paths
    // Iterate through all paths and write to file
    size_t path_count = 0;
    // gbz.graph.for_each_path_handle([&](const path_handle_t& path) {
    //     std::string path_name = gbz.graph.get_path_name(path);
    //     std::cerr << "Path name ( " << path_count << "): " << path_name << std::endl;
    //     output << path_name << std::endl;
    //     path_count++;
    // });


    if (gbz.index.metadata.hasPathNames()) {
        auto gbwt_reference_samples = gbwtgraph::parse_reference_samples_tag(gbz.index);
        for (size_t i = 0; i < gbz.index.metadata.paths(); i++) {
            PathSense sense = gbwtgraph::get_path_sense(gbz.index, i, gbwt_reference_samples);
            std::string path_name = gbwtgraph::compose_path_name(gbz.index, i, sense);
            std::cout << path_name << std::endl;
            output << path_name << std::endl;
            path_count++;
        }
    }


    output.close();

    std::cerr << "Extracted " << path_count << " paths to: " << output_file << std::endl;

    return 0;
}