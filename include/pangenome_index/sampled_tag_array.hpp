#ifndef PANGENOME_INDEX_SAMPLED_TAG_ARRAY_HPP
#define PANGENOME_INDEX_SAMPLED_TAG_ARRAY_HPP

#include <cstdint>
#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_gmr.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/util.hpp>
#include <mutex>
#include <gbwtgraph/utils.h>

namespace panindexer {

    class SampledTagArray {
    public:
        SampledTagArray() = default;

        // Build from a stream of runs: for each input run (pos_t, length),
        // emit value = encode(node_id,is_rev) if offset==0, else GAP_CODE (0), merging consecutive runs with equal value.
        void build_from_runs(const std::vector<std::pair<handlegraph::pos_t, uint64_t>>& runs, size_t bwt_size);

        // Build from a callback enumerator (yields many runs without materializing all of them)
        void build_from_enumerator(const std::function<void(const std::function<void(handlegraph::pos_t,uint64_t)>&)>& enumerator,
                                   size_t bwt_size);

        // Serialization
        void serialize(std::ostream& out) const;
        void load(std::istream& in);

        // Accessors
        inline const sdsl::wt_gmr<>& values() const { return sampled_values; }
        inline const sdsl::sd_vector<>& run_starts() const { return bwt_intervals; }

        // Lazy support initialization for run_starts()
        inline void ensure_run_rank() const {
            std::call_once(run_rank_once, [&]() {
                sdsl::util::init_support(run_rank_support, &bwt_intervals);
            });
        }

        inline void ensure_run_select() const {
            std::call_once(run_select_once, [&]() {
                sdsl::util::init_support(run_select_support, &bwt_intervals);
            });
        }

        // Helpers for queries
        inline size_t total_runs() const { return sampled_values.size(); }

        // Return run id that contains BWT position pos
        inline size_t run_id_at(size_t pos) const {
            ensure_run_rank();
            return run_rank_support(pos);
        }

        // Return [start,end] BWT span for run_id
        inline std::pair<size_t,size_t> run_span(size_t run_id) const {
            ensure_run_select();
            size_t start = run_select_support(run_id + 1);
            size_t end;
            if (run_id + 1 < total_runs()) {
                end = run_select_support(run_id + 2) - 1;
            } else {
                end = bwt_intervals.size() - 1;
            }
            return { start, end };
        }

        // Return encoded tag value for a run
        inline uint64_t run_value(size_t run_id) const { return sampled_values[run_id]; }

        // Encode (node_id,is_rev) into integer code; 0 reserved for gaps
        static inline uint64_t encode_value(int64_t node_id, bool is_rev) {
            // node_id >= 1; map to 1.. via shift and set bit0 for strand
            return 1 + ((static_cast<uint64_t>(node_id - 1) << 1) | static_cast<uint64_t>(is_rev));
        }

    private:
        sdsl::wt_gmr<> sampled_values; // one value per run (0 for gap, otherwise encode(node_id,is_rev))
        sdsl::sd_vector<> bwt_intervals; // 1 at BWT positions where a run starts

        // Lazy supports for run_starts
        mutable sdsl::sd_vector<>::rank_1_type run_rank_support;
        mutable sdsl::sd_vector<>::select_1_type run_select_support;
        mutable std::once_flag run_rank_once;
        mutable std::once_flag run_select_once;
    };

}

#endif // PANGENOME_INDEX_SAMPLED_TAG_ARRAY_HPP


