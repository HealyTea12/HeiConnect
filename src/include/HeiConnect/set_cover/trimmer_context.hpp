#pragma once

#include "HeiConnect/set_cover/common.hpp"

#include <memory>
#include <unordered_set>
#include <algorithm>

// Bit-packed trimmer context: specialized for bit / pseudo representations
// that provide `get_col(set_index, col_index)` and `get_n_cols()` helpers.
template <typename SetCoverType, typename SolutionType>
    requires ForEachElementCon<SetCoverType>
class BitPackedTrimmerContext
{
public:
    explicit BitPackedTrimmerContext(std::shared_ptr<const SetCoverType> set_cover)
        : m_set_cover(std::move(set_cover))
    {
    }

    // Doesn't need to keep any state, but must satisfy concept
    void add_set(size_t set_index)
    {
    }

    void remove_set(size_t set_index)
    {
    }

    bool can_remove(size_t set_index, const SolutionType &solution) const
    {
        // For bit-packed representations we can inspect columns
        // and test whether any bits remain after removing other sets.
        const size_t n_cols = m_set_cover->get_n_cols();
        for (size_t k = 0; k < n_cols; ++k)
        {
            ull set_coverage = m_set_cover->get_col(set_index, k);
            for (const auto &other_set_index : solution.get_solution())
            {
                if (other_set_index == set_index)
                    continue;
                set_coverage &= ~(m_set_cover->get_col(other_set_index, k));
                if (set_coverage == 0)
                    break;
            }
            if (set_coverage != 0)
            {
                return false;
            }
        }
        return true;
    }

private:
    std::shared_ptr<const SetCoverType> m_set_cover;
};

// Keep a placeholder specialization for SetCoverDouble -> callers can
// implement trimming for the double-representation if needed.
// Need to think how to implement this best
/*
template <class link_node_T, class link_edge_T, class link_weight_T>
class BitPackedTrimmerContext<
    SetCoverDouble<
        SetCoverPseudo<link_node_T, link_edge_T>,
        SetCoverCyc<link_node_T, link_edge_T, link_weight_T>>, USSolution>
{
public:
    using SetCoverT = SetCoverDouble<
        SetCoverPseudo<link_node_T, link_edge_T>,
        SetCoverCyc<link_node_T, link_edge_T, link_weight_T>>;

    explicit BitPackedTrimmerContext(std::shared_ptr<const SetCoverT>) {}

    bool can_remove(typename SetCoverT::SetID, const USSolution &)
    {
        throw std::logic_error("Trimming not implemented for SetCoverDouble with SetCoverCyc");
    }

    void add_set(typename SetCoverT::SetID) {}
    void remove_set(typename SetCoverT::SetID) {}
    size_t get_total_covered_elements() const { return 0; }
    void reset() {}
};
*/
