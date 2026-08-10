#include "pangenome_index/tag_head_samples.hpp"

#include <sdsl/io.hpp>
#include <sdsl/util.hpp>

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

namespace panindexer {

// Out-of-class definitions for the constexpr static members (ODR-use safety).
constexpr char     TagHeadSamples::MAGIC[8];
constexpr uint64_t TagHeadSamples::VERSION;
constexpr uint64_t TagHeadSamples::UNKNOWN;

TagHeadSamples::TagHeadSamples() = default;

void TagHeadSamples::rebind_supports() {
    sdsl::util::init_support(kept_rank_, &kept_marker_);
}

void TagHeadSamples::load(std::istream& in) {
    if (!in) throw std::runtime_error("TagHeadSamples::load: bad istream");

    char     m[8];
    uint64_t ver = 0, s_u = 0, ntr = 0, kc = 0, nb = 0;
    in.read(m, sizeof(m));
    in.read(reinterpret_cast<char*>(&ver), sizeof(ver));
    in.read(reinterpret_cast<char*>(&s_u), sizeof(s_u));
    in.read(reinterpret_cast<char*>(&ntr), sizeof(ntr));
    in.read(reinterpret_cast<char*>(&kc),  sizeof(kc));
    in.read(reinterpret_cast<char*>(&nb),  sizeof(nb));
    if (!in) throw std::runtime_error("TagHeadSamples::load: header read failed");

    if (std::memcmp(m, MAGIC, sizeof(MAGIC)) != 0) {
        throw std::runtime_error("TagHeadSamples::load: bad magic (not a .tag_samples file)");
    }
    if (ver != VERSION) {
        throw std::runtime_error("TagHeadSamples::load: unsupported version "
                                 + std::to_string(ver) + " (expected "
                                 + std::to_string(VERSION) + ")");
    }
    if (s_u == 0) {
        throw std::runtime_error("TagHeadSamples::load: invalid s = 0");
    }

    s_            = static_cast<std::size_t>(s_u);
    num_tag_runs_ = static_cast<std::size_t>(ntr);
    kept_count_   = static_cast<std::size_t>(kc);
    bwt_size_     = static_cast<std::size_t>(nb);

    kept_marker_.load(in);
    sa_values_.load(in);
    if (!in) throw std::runtime_error("TagHeadSamples::load: body read failed");

    // Structural consistency checks.
    if (kept_marker_.size() != num_tag_runs_) {
        throw std::runtime_error("TagHeadSamples::load: kept_marker size ("
                                 + std::to_string(kept_marker_.size())
                                 + ") != num_tag_runs (" + std::to_string(num_tag_runs_) + ")");
    }
    if (sa_values_.size() != kept_count_) {
        throw std::runtime_error("TagHeadSamples::load: sa_values length ("
                                 + std::to_string(sa_values_.size())
                                 + ") != kept_count (" + std::to_string(kept_count_) + ")");
    }

    rebind_supports();

    // Verify kept_marker 1-count matches kept_count (belt-and-suspenders).
    std::size_t actual_ones = (kept_marker_.size() == 0)
                                  ? 0
                                  : kept_rank_(kept_marker_.size());
    if (actual_ones != kept_count_) {
        throw std::runtime_error("TagHeadSamples::load: kept_marker 1-count ("
                                 + std::to_string(actual_ones)
                                 + ") != kept_count (" + std::to_string(kept_count_) + ")");
    }
}

std::size_t TagHeadSamples::memory_bytes() const {
    std::size_t total = 0;
    total += sdsl::size_in_bytes(kept_marker_);
    total += sdsl::size_in_bytes(kept_rank_);
    total += sdsl::size_in_bytes(sa_values_);
    return total;
}

} // namespace panindexer
