/**
 * Translation Table 1 implementation.
 * Converts (named path, global interval) → list of (path_id, local interval).
 */

#include "pangenome_index/translation_tables.hpp"
#include <algorithm>
#include <cassert>
#include <istream>
#include <ostream>
#include <stdexcept>

namespace panindexer {

void TranslationTable1::add_subpath(const std::string& path_name,
                                     size_t path_id,
                                     size_t subpath_start,
                                     size_t length) {
    if (name_to_subpaths_.find(path_name) == name_to_subpaths_.end()) {
        names_.push_back(path_name);
    }
    name_to_subpaths_[path_name].push_back(
        SubpathInfo{path_id, subpath_start, length});
}

std::vector<PathInterval> TranslationTable1::lookup(const std::string& path_name,
                                                    size_t global_start,
                                                    size_t global_end) const {
    auto it = name_to_subpaths_.find(path_name);
    if (it == name_to_subpaths_.end()) {
        return {};
    }
    const std::vector<SubpathInfo>& subpaths = it->second;
    std::vector<PathInterval> result;
    for (const SubpathInfo& sp : subpaths) {
        size_t subpath_end = sp.subpath_start + sp.length;
        // Overlap [global_start, global_end) with [sp.subpath_start, subpath_end)
        size_t overlap_start = std::max(global_start, sp.subpath_start);
        size_t overlap_end = std::min(global_end, subpath_end);
        if (overlap_start >= overlap_end) {
            continue;
        }
        PathInterval pi;
        pi.path_id = sp.path_id;
        pi.start = overlap_start - sp.subpath_start;
        pi.end = overlap_end - sp.subpath_start;
        result.push_back(pi);
    }
    return result;
}

std::vector<PathInterval> TranslationTable1::lookup_by_name_id(size_t name_id,
                                                               size_t global_start,
                                                               size_t global_end) const {
    if (name_id >= names_.size()) {
        return {};
    }
    return lookup(names_[name_id], global_start, global_end);
}

std::vector<std::string> TranslationTable1::names() const {
    return names_;
}

std::vector<SubpathInfo> TranslationTable1::subpaths(const std::string& path_name) const {
    auto it = name_to_subpaths_.find(path_name);
    if (it == name_to_subpaths_.end()) {
        return {};
    }
    return it->second;
}

namespace {

const uint32_t TABLE1_MAGIC = 0x54543100;  // "TT1\0"
const uint32_t TABLE1_VERSION = 1;

void write_uint32(std::ostream& out, uint32_t x) {
    out.put(static_cast<char>(x & 0xff));
    out.put(static_cast<char>((x >> 8) & 0xff));
    out.put(static_cast<char>((x >> 16) & 0xff));
    out.put(static_cast<char>((x >> 24) & 0xff));
}

void write_uint64(std::ostream& out, uint64_t x) {
    for (int i = 0; i < 8; ++i) {
        out.put(static_cast<char>(x & 0xff));
        x >>= 8;
    }
}

uint32_t read_uint32(std::istream& in) {
    uint32_t x = 0;
    for (int i = 0; i < 4; ++i) {
        int c = in.get();
        if (c == std::char_traits<char>::eof()) {
            throw std::runtime_error("TranslationTable1::load: unexpected EOF");
        }
        x |= static_cast<uint32_t>(static_cast<unsigned char>(c)) << (i * 8);
    }
    return x;
}

uint64_t read_uint64(std::istream& in) {
    uint64_t x = 0;
    for (int i = 0; i < 8; ++i) {
        int c = in.get();
        if (c == std::char_traits<char>::eof()) {
            throw std::runtime_error("TranslationTable1::load: unexpected EOF");
        }
        x |= static_cast<uint64_t>(static_cast<unsigned char>(c)) << (i * 8);
    }
    return x;
}

void write_string(std::ostream& out, const std::string& s) {
    write_uint64(out, s.size());
    out.write(s.data(), static_cast<std::streamsize>(s.size()));
}

std::string read_string(std::istream& in) {
    uint64_t len = read_uint64(in);
    std::string s(len, '\0');
    in.read(&s[0], static_cast<std::streamsize>(len));
    if (!in) {
        throw std::runtime_error("TranslationTable1::load: failed to read string");
    }
    return s;
}

} // namespace

void TranslationTable1::serialize(std::ostream& out) const {
    write_uint32(out, TABLE1_MAGIC);
    write_uint32(out, TABLE1_VERSION);
    write_uint64(out, names_.size());
    for (const std::string& name : names_) {
        write_string(out, name);
        const std::vector<SubpathInfo>& subpaths = name_to_subpaths_.at(name);
        write_uint64(out, subpaths.size());
        for (const SubpathInfo& sp : subpaths) {
            write_uint64(out, sp.path_id);
            write_uint64(out, sp.subpath_start);
            write_uint64(out, sp.length);
        }
    }
}

void TranslationTable1::load(std::istream& in) {
    uint32_t magic = read_uint32(in);
    if (magic != TABLE1_MAGIC) {
        throw std::runtime_error("TranslationTable1::load: invalid magic (not a Table1 file?)");
    }
    uint32_t version = read_uint32(in);
    if (version != TABLE1_VERSION) {
        throw std::runtime_error("TranslationTable1::load: unsupported version");
    }
    name_to_subpaths_.clear();
    names_.clear();
    uint64_t num_names = read_uint64(in);
    for (uint64_t i = 0; i < num_names; ++i) {
        std::string name = read_string(in);
        uint64_t num_subpaths = read_uint64(in);
        std::vector<SubpathInfo> subpaths;
        subpaths.reserve(num_subpaths);
        for (uint64_t j = 0; j < num_subpaths; ++j) {
            SubpathInfo sp;
            sp.path_id = read_uint64(in);
            sp.subpath_start = read_uint64(in);
            sp.length = read_uint64(in);
            subpaths.push_back(sp);
        }
        names_.push_back(name);
        name_to_subpaths_[name] = std::move(subpaths);
    }
}

// ============================================================
// TranslationTable2 implementation
// ============================================================

void TranslationTable2::add_mapping(size_t src_path_id,
                                     const std::string& tgt_haplotype,
                                     const IntervalMapping& mapping) {
    Key k{src_path_id, tgt_haplotype};
    entries_[k].push_back(mapping);
}

void TranslationTable2::finalize() {
    for (auto& kv : entries_) {
        auto& segs = kv.second;
        std::sort(segs.begin(), segs.end(),
                  [](const IntervalMapping& a, const IntervalMapping& b) {
                      return a.src_start < b.src_start;
                  });
        // Overlapping source intervals are allowed: src[a,b] and src[c,d] can both
        // map to different tgt_path_ids (e.g. src[0,100]->i, src[10,90]->j).
    }
}

std::vector<TargetInterval> TranslationTable2::lookup(size_t src_path_id,
                                                       const std::string& tgt_haplotype,
                                                       size_t local_start,
                                                       size_t local_end) const {
    Key k{src_path_id, tgt_haplotype};
    auto it = entries_.find(k);
    if (it == entries_.end()) {
        return {};
    }
    const auto& segs = it->second;
    if (segs.empty() || local_start >= local_end) {
        return {};
    }

    // Binary search for first segment whose src_end > local_start
    size_t lo = 0, hi = segs.size();
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        if (segs[mid].src_end <= local_start) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    std::vector<TargetInterval> result;
    for (size_t i = lo; i < segs.size(); ++i) {
        const IntervalMapping& seg = segs[i];
        if (seg.src_start >= local_end) break;  // past query range

        size_t clip_src_start = std::max(local_start, seg.src_start);
        size_t clip_src_end   = std::min(local_end,   seg.src_end);
        if (clip_src_start >= clip_src_end) continue;

        TargetInterval ti;
        ti.tgt_path_id = seg.tgt_path_id;
        result.push_back(ti);
    }
    return result;
}

size_t TranslationTable2::total_segments() const {
    size_t total = 0;
    for (const auto& kv : entries_) {
        total += kv.second.size();
    }
    return total;
}

std::vector<std::pair<size_t, std::string>> TranslationTable2::keys() const {
    std::vector<std::pair<size_t, std::string>> result;
    result.reserve(entries_.size());
    for (const auto& kv : entries_) {
        result.emplace_back(kv.first.src_path_id, kv.first.tgt_haplotype);
    }
    return result;
}

std::vector<IntervalMapping> TranslationTable2::segments(size_t src_path_id,
                                                          const std::string& tgt_haplotype) const {
    Key k{src_path_id, tgt_haplotype};
    auto it = entries_.find(k);
    if (it == entries_.end()) return {};
    return it->second;
}

namespace {

const uint32_t TABLE2_MAGIC   = 0x54543200;  // "TT2\0"
const uint32_t TABLE2_VERSION = 2;  // v2: no tgt_start/tgt_end in IntervalMapping

}  // namespace

void TranslationTable2::serialize(std::ostream& out) const {
    write_uint32(out, TABLE2_MAGIC);
    write_uint32(out, TABLE2_VERSION);
    write_uint64(out, static_cast<uint64_t>(entries_.size()));
    for (const auto& kv : entries_) {
        write_uint64(out, static_cast<uint64_t>(kv.first.src_path_id));
        write_string(out, kv.first.tgt_haplotype);
        const auto& segs = kv.second;
        write_uint64(out, static_cast<uint64_t>(segs.size()));
        for (const IntervalMapping& seg : segs) {
            write_uint64(out, static_cast<uint64_t>(seg.src_start));
            write_uint64(out, static_cast<uint64_t>(seg.src_end));
            write_uint64(out, static_cast<uint64_t>(seg.tgt_path_id));
        }
    }
}

void TranslationTable2::load(std::istream& in) {
    uint32_t magic = read_uint32(in);
    if (magic != TABLE2_MAGIC) {
        throw std::runtime_error("TranslationTable2::load: invalid magic (not a Table2 file?)");
    }
    uint32_t version = read_uint32(in);
    if (version != 1 && version != TABLE2_VERSION) {
        throw std::runtime_error("TranslationTable2::load: unsupported version");
    }
    entries_.clear();
    uint64_t num_keys = read_uint64(in);
    for (uint64_t i = 0; i < num_keys; ++i) {
        Key k;
        k.src_path_id  = static_cast<size_t>(read_uint64(in));
        k.tgt_haplotype = read_string(in);
        uint64_t num_segs = read_uint64(in);
        std::vector<IntervalMapping> segs;
        segs.reserve(static_cast<size_t>(num_segs));
        for (uint64_t j = 0; j < num_segs; ++j) {
            IntervalMapping seg;
            seg.src_start   = static_cast<size_t>(read_uint64(in));
            seg.src_end     = static_cast<size_t>(read_uint64(in));
            seg.tgt_path_id = static_cast<size_t>(read_uint64(in));
            if (version == 1) {
                (void)read_uint64(in);  // discard tgt_start
                (void)read_uint64(in);  // discard tgt_end
            }
            segs.push_back(seg);
        }
        entries_[k] = std::move(segs);
    }
}

} // namespace panindexer
