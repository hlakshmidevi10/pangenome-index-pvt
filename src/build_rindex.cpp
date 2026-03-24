//
// Created by seeskand on 12/9/24.
//

#include "pangenome_index/r-index.hpp"
#include <iostream>
#include <cstring>

using namespace panindexer;

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <rlbwt_file> [> output.ri]\n"
              << "\n"
              << "Build the r-index (FastLocate) from an RLBWT file and write it to stdout.\n"
              << "Redirect stdout to save the index (e.g. " << prog << " text.rl_bwt > rindex.ri).\n"
              << "\n"
              << "Arguments:\n"
              << "  <rlbwt_file>    Input RLBWT file (e.g. .rl_bwt)\n"
              << "\n"
              << "Options:\n"
              << "  -h, --help      Print this usage and exit\n"
              << std::endl;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }
    if (std::strcmp(argv[1], "-h") == 0 || std::strcmp(argv[1], "--help") == 0) {
        usage(argv[0]);
        return 0;
    }

    std::string rlbwt_file = argv[1];

    FastLocate idx(rlbwt_file);
    idx.serialize_encoded(std::cout);

    return 0;
}