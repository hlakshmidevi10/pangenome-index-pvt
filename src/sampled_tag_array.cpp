#include "pangenome_index/sampled_tag_array.hpp"
#include "pangenome_index/tag_arrays.hpp"
#include <sdsl/util.hpp>

namespace panindexer {

    using handlegraph::pos_t;

    static inline uint64_t encode_val_from_pos(pos_t p) {
        if (gbwtgraph::offset(p) != 0) return 0; // gap within node
        if (gbwtgraph::id(p) == 0) return 0;     // treat special zero-tag as gap
        return SampledTagArray::encode_value(gbwtgraph::id(p), gbwtgraph::is_rev(p));
    }

    // Helper function to construct wt_gmr from values
    // wt_gmr uses grammar-based compression which is more memory-efficient and stable on macOS
    static void construct_wt_gmr_from_values(sdsl::wt_gmr<sdsl::int_vector<>, sdsl::inv_multi_perm_support<4, sdsl::int_vector<>>>& target, const std::vector<uint64_t>& values) {
        if (values.empty()) {
            target = sdsl::wt_gmr<sdsl::int_vector<>, sdsl::inv_multi_perm_support<4, sdsl::int_vector<>>>();
            return;
        }
        
        uint64_t maxv = 0;
        for (auto v : values) { if (v > maxv) maxv = v; }
        size_t width = std::max<size_t>(1, sdsl::bits::hi(maxv) + 1);
        
        // Create int_vector and populate it
        sdsl::int_vector<> iv(values.size(), 0, width);
        for (size_t i = 0; i < values.size(); ++i) {
            iv[i] = values[i];
        }
        
        // Construct wt_gmr directly - it should handle memory more reliably than wm_int
        sdsl::construct_im(target, iv);
        
        // iv will be destroyed when it goes out of scope
        // target should have its own independent copy of the data
    }

    void SampledTagArray::build_from_runs(const std::vector<std::pair<pos_t, uint64_t>>& runs, size_t bwt_size) {
        // First, transform runs into (start position markers) using cumulative lengths, and values merged.
        // We'll collect run_starts positions and run values in order.
        std::vector<uint64_t> values;
        values.reserve(runs.size());

        std::vector<uint64_t> starts;
        starts.reserve(runs.size());

        uint64_t cum = 0;
        uint64_t last_val = std::numeric_limits<uint64_t>::max();

        for (const auto& pr : runs) {
            pos_t p = pr.first;
            uint64_t len = pr.second;
            uint64_t val = encode_val_from_pos(p);
            if (values.empty() || val != last_val) {
                starts.push_back(cum);
                values.push_back(val);
                last_val = val;
            }
            cum += len;
        }

        // Determine if first run is gap (false if gap, true if normal tag)
        first_run_is_gap = (values.empty() || values[0] == 0);

        // Build bwt_intervals with zero-length gap runs inserted between adjacent normal tags
        // Only store non-gap values in wt_gmr
        std::vector<uint64_t> bwt_starts; // All run starts including gaps
        std::vector<uint64_t> non_gap_values; // Only non-gap values for wt_gmr
        
        for (size_t i = 0; i < values.size(); ++i) {
            uint64_t val = values[i];
            uint64_t start = starts[i];
            
            // Add start position for this run
            bwt_starts.push_back(start);
            
            // If this is a non-gap value, store it in wt_gmr
            if (val != 0) {
                non_gap_values.push_back(val);
                
                // Insert zero-length gap run between this and next run if next run is also non-gap
                if (i + 1 < values.size() && values[i + 1] != 0) {
                    // Insert gap run at the start position of next run (zero-length gap)
                    bwt_starts.push_back(starts[i + 1]);
                }
            }
        }

        // Build sd_vector for run starts with multiset=true
        {
            sdsl::sd_vector_builder builder(bwt_size + 1, bwt_starts.size(), true); // multiset=true
            for (uint64_t s : bwt_starts) builder.set(s);
            bwt_intervals = sdsl::sd_vector<>(builder);
        }

        // Build wt_gmr from only non-gap values
        construct_wt_gmr_from_values(sampled_values, non_gap_values);
    }

    void SampledTagArray::build_from_enumerator(const std::function<void(const std::function<void(pos_t,uint64_t)>&)>& enumerator,
                                                size_t bwt_size) {
        std::vector<uint64_t> values;
        std::vector<uint64_t> starts;
        uint64_t cum = 0;
        uint64_t last_val = std::numeric_limits<uint64_t>::max();

        auto sink = [&](pos_t p, uint64_t len) {
            uint64_t val = encode_val_from_pos(p);
            if (values.empty() || val != last_val) {
                starts.push_back(cum);
                values.push_back(val);
                last_val = val;
            }
            cum += len;
        };

        enumerator(sink);
        
        // Determine if first run is gap (with bounds check for empty input)
        first_run_is_gap = (values.empty() || values[0] == 0);

        // Build bwt_intervals with zero-length gap runs inserted between adjacent normal tags
        // Only store non-gap values in wt_gmr
        std::vector<uint64_t> bwt_starts; // All run starts including gaps
        std::vector<uint64_t> non_gap_values; // Only non-gap values for wt_gmr
        
        for (size_t i = 0; i < values.size(); ++i) {
            uint64_t val = values[i];
            uint64_t start = starts[i];
            
            // Add start position for this run
            bwt_starts.push_back(start);
            
            // If this is a non-gap value, store it in wt_gmr
            if (val != 0) {
                non_gap_values.push_back(val);
                
                // Insert zero-length gap run between this and next run if next run is also non-gap
                if (i + 1 < values.size() && values[i + 1] != 0) {
                    // Insert gap run at the start position of next run (zero-length gap)
                    bwt_starts.push_back(starts[i + 1]);
                }
            }
        }
        
        // Build sd_vector for run starts with multiset=true
        {
            sdsl::sd_vector_builder builder(bwt_size + 1, bwt_starts.size(), true); // multiset=true
            for (uint64_t s : bwt_starts) builder.set(s);
            bwt_intervals = sdsl::sd_vector<>(builder);
        }
        
        // Build wt_gmr from only non-gap values
        construct_wt_gmr_from_values(sampled_values, non_gap_values);
        
        std::cerr << "Finished building sampled_tag_array" << std::endl;
    }

    void SampledTagArray::serialize(std::ostream& out) const {
        // On macOS, ensure supports are not accessed during serialization
        // Serialize the structures directly without accessing lazy supports
        sdsl::serialize(sampled_values, out);
        sdsl::serialize(bwt_intervals, out);
        sdsl::write_member(first_run_is_gap, out);
    }

    void SampledTagArray::load(std::istream& in) {
        sampled_values.load(in);
        bwt_intervals.load(in);
        sdsl::read_member(first_run_is_gap, in);
    }

}


