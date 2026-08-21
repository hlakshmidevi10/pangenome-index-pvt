// probe_blocks: read-only diagnostic for the block-count off-by-one.
//
// Loads a .ri and reports whether `blocks` / `blocks_encoded_start_bits`
// (sized floor(total_runs/block_size)+1 by the constructor) disagrees with
// the number of 1-bits in `blocks_start_pos` (sized ceil(total_runs/block_size)
// by the sd_vector_builder). They diverge exactly when
// total_runs % block_size == 0, leaving a phantom trailing block.
//
// Temporary; delete once the off-by-one is fixed.

#include "pangenome_index/r-index.hpp"

#include <sdsl/sd_vector.hpp>

#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

using namespace panindexer;

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "usage: probe_blocks <ri_file>\n"; return 2; }

    FastLocate ri;
    {
        std::ifstream in(argv[1], std::ios::binary);
        if (!in) { std::cerr << "ERROR: cannot open " << argv[1] << "\n"; return 1; }
        ri.load_encoded(in);
    }

    const size_t total_runs = ri.tot_runs();
    const size_t nb         = ri.num_blocks();
    const size_t bs         = ri.encoded_block_size ? ri.encoded_block_size : ri.block_size;

    sdsl::sd_vector<>::rank_1_type rk(&ri.blocks_start_pos);
    const size_t ones = rk(ri.blocks_start_pos.size());

    std::cout << "ri                        = " << argv[1] << "\n"
              << "tot_runs()                = " << total_runs << "\n"
              << "samples.size()            = " << ri.samples.size() << "\n"
              << "samples.width()           = " << static_cast<unsigned>(ri.samples.width()) << " bits\n"
              << "effective block_size      = " << bs << "\n"
              << "total_runs % block_size   = " << (total_runs % bs) << "\n"
              << "C.size()                  = " << ri.C.size() << "\n"
              << "--- block accounting ---\n"
              << "num_blocks()              = " << nb << "\n"
              << "  floor(runs/bs) + 1      = " << (total_runs / bs + 1)
              << "   <-- constructor r-index.cpp:1127\n"
              << "  ceil(runs/bs)           = " << (total_runs / bs + (total_runs % bs != 0))
              << "   <-- sd_builder r-index.cpp:1143\n"
              << "blocks_start_pos ones     = " << ones << "\n"
              << "MISMATCH (num_blocks-ones)= " << (static_cast<long long>(nb) - static_cast<long long>(ones))
              << "\n";

    std::cout << "--- encoded stream tail ---\n"
              << "blocks_encoded_stream size= " << ri.blocks_encoded_stream.size() << "\n";
    for (size_t b = (nb >= 3 ? nb - 3 : 0); b < nb; ++b) {
        std::cout << "  start_bits[" << b << "] = " << ri.blocks_encoded_start_bits[b] << "\n";
    }

    std::cout << "--- decoding trailing blocks ---\n";
    for (size_t b = (nb >= 3 ? nb - 3 : 0); b < nb; ++b) {
        std::vector<std::pair<size_t, size_t>> runs;
        ri.get_block_runs(b, runs);
        std::cout << "  get_block_runs(" << b << ") -> " << runs.size() << " runs:";
        for (size_t i = 0; i < runs.size() && i < 12; ++i) {
            std::cout << " (sym=" << runs[i].first << ",len=" << runs[i].second << ")";
        }
        std::cout << "\n";
    }

    std::cout << "--- verdict ---\n";
    if (nb > ones) {
        std::cout << "  blocks_start_select_1(" << nb << ") would read low["
                  << (nb - 1) << "] of size " << ones << "  => OUT OF BOUNDS => ABORT\n";
    } else {
        std::cout << "  block accounting consistent; no phantom block.\n";
    }
    return 0;
}
