#include "pangenome_index/light_tag_index.hpp"

#include "pangenome_index/tag_arrays.hpp"

#include <sdsl/io.hpp>
#include <sdsl/util.hpp>

#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace panindexer {

// MAGIC[8] is defined inline in the header (constexpr char array). We also need
// an out-of-class definition for the C++17 ODR-use case; this is a no-op for
// inline constexpr but kept for safety with non-inline interpretations.
constexpr char     LightTagIndex::MAGIC[8];
constexpr uint64_t LightTagIndex::VERSION;

LightTagIndex::LightTagIndex() = default;

void LightTagIndex::rebind_supports() {
    // sdsl rank/select supports hold a pointer back to the bit vector. After
    // load() or set_bwt_intervals() we must rebind to the current member.
    sdsl::util::init_support(bwt_intervals_rank,  &bwt_intervals);
    sdsl::util::init_support(bwt_intervals_select, &bwt_intervals);
    // Total 1-count = total number of tag runs.
    num_runs_cached = (bwt_intervals.size() == 0)
                          ? 0
                          : bwt_intervals_rank(bwt_intervals.size());
}

void LightTagIndex::set_bwt_intervals(sdsl::sd_vector<>&& v) {
    bwt_intervals = std::move(v);
    rebind_supports();
}

void LightTagIndex::save(std::ostream& out) const {
    if (!out) throw std::runtime_error("LightTagIndex::save: bad ostream");
    // Header: magic + version + bwt_size.
    out.write(MAGIC, sizeof(MAGIC));
    uint64_t ver = VERSION;
    out.write(reinterpret_cast<const char*>(&ver), sizeof(ver));
    uint64_t bsz = static_cast<uint64_t>(bwt_intervals.size());
    out.write(reinterpret_cast<const char*>(&bsz), sizeof(bsz));
    if (!out) throw std::runtime_error("LightTagIndex::save: write failed (header)");
    // Body: sd_vector serialization.
    bwt_intervals.serialize(out);
    if (!out) throw std::runtime_error("LightTagIndex::save: write failed (body)");
}

void LightTagIndex::load(std::istream& in) {
    if (!in) throw std::runtime_error("LightTagIndex::load: bad istream");
    char     m[8];
    uint64_t ver = 0;
    uint64_t bsz = 0;
    in.read(m, sizeof(m));
    in.read(reinterpret_cast<char*>(&ver), sizeof(ver));
    in.read(reinterpret_cast<char*>(&bsz), sizeof(bsz));
    if (!in) throw std::runtime_error("LightTagIndex::load: header read failed");
    if (std::memcmp(m, MAGIC, sizeof(MAGIC)) != 0) {
        throw std::runtime_error("LightTagIndex::load: bad magic (not a .ltags file)");
    }
    if (ver != VERSION) {
        throw std::runtime_error("LightTagIndex::load: unsupported version "
                                 + std::to_string(ver) + " (expected "
                                 + std::to_string(VERSION) + ")");
    }
    bwt_intervals.load(in);
    if (!in) throw std::runtime_error("LightTagIndex::load: sd_vector load failed");
    if (bwt_intervals.size() != bsz) {
        throw std::runtime_error("LightTagIndex::load: header bwt_size ("
                                 + std::to_string(bsz)
                                 + ") != sd_vector size ("
                                 + std::to_string(bwt_intervals.size()) + ")");
    }
    rebind_supports();
}

void LightTagIndex::build_from_compact_tags_file(const std::string& compact_tags_path) {
    std::ifstream in(compact_tags_path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("LightTagIndex::build_from_compact_tags_file: "
                                 "cannot open " + compact_tags_path);
    }
    TagArray tags;
    tags.load_compressed_tags_compact(in);
    // Copy bwt_intervals out (sd_vector is copy-constructible).
    sdsl::sd_vector<> copy = tags.get_bwt_intervals();
    set_bwt_intervals(std::move(copy));
}

std::size_t LightTagIndex::memory_bytes() const {
    std::size_t total = 0;
    total += sdsl::size_in_bytes(bwt_intervals);
    total += sdsl::size_in_bytes(bwt_intervals_rank);
    total += sdsl::size_in_bytes(bwt_intervals_select);
    return total;
}

} // namespace panindexer
