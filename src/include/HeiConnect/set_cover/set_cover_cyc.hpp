#pragma once

#include "HeiConnect/set_cover/common.hpp"

#include <ranges>

template <typename node_T, typename cycle_id_T>
struct CycleCross
{
    cycle_id_T cycle_id;
    node_T v1;
    node_T v2;
};

template <class link_node_T, class link_edge_T, class link_weight_T>
    requires std::integral<link_node_T> && std::integral<link_edge_T>
class SetCoverCyc
{
public:
    using cycle_pos_T = int;
    using cycle_id_T = int;
    using ull = unsigned long long;
    using ElementID = size_t;
    using SetID = size_t;
    using SetCost = double;
    struct Link
    {
        link_node_T u;
        link_node_T v;
        SetCost weight;
    };

    class SetElementsIterator
    {
    public:
        using value_type = size_t;
        using difference_type = std::ptrdiff_t;
        using iterator_concept = std::input_iterator_tag;
        using iterator_category = std::input_iterator_tag;

        SetElementsIterator(const SetCoverCyc &owner, size_t set_index, bool is_end)
            : m_owner(owner), m_set_index(set_index), m_is_end(is_end)
        {
            if (m_is_end)
            {
                return;
            }

            const auto &cyc = m_owner.get();
            if (m_set_index >= cyc.links.size())
            {
                m_is_end = true;
                return;
            }

            advance();
        }

        SetElementsIterator &operator++()
        {
            advance();
            return *this;
        }

        SetElementsIterator operator++(int)
        {
            SetElementsIterator tmp = *this;
            ++(*this);
            return tmp;
        }

        size_t operator*() const
        {
            return m_current_element;
        }

        bool operator!=(const SetElementsIterator &other) const
        {
            return m_is_end != other.m_is_end;
        }

        friend bool operator==(const SetElementsIterator &a, const SetElementsIterator &b)
        {
            return !(a != b);
        }

    private:
        enum class Phase
        {
            Cycle,
            End
        };

        bool prepare_cycle_cross()
        {
            const auto &cyc = m_owner.get();
            const auto &crosses = cyc.m_cycle_crosses[m_set_index];
            while (m_cross_idx < crosses.size())
            {
                const auto &[cid, a, b] = crosses[m_cross_idx];
                const size_t cycle_idx = static_cast<size_t>(cid);
                if (cycle_idx >= cyc.cycle_sizes.size())
                {
                    ++m_cross_idx;
                    continue;
                }

                m_cycle_size = cyc.cycle_sizes[cycle_idx];
                m_a = static_cast<size_t>(a);
                m_b = static_cast<size_t>(b);
                m_cycle_offset = cyc.m_cycle_element_offsets[cycle_idx];

                m_in_second_region = false;
                m_j = 0;
                m_k = m_a;

                if (m_a > m_cycle_size)
                {
                    m_a = m_cycle_size;
                }
                if (m_b > m_cycle_size)
                {
                    m_b = m_cycle_size;
                }
                if (m_b < m_a)
                {
                    std::swap(m_a, m_b);
                }
                return true;
            }
            return false;
        }

        bool advance_cycle()
        {
            const auto &cyc = m_owner.get();
            const auto &crosses = cyc.m_cycle_crosses[m_set_index];

            while (true)
            {
                if (m_cross_idx >= crosses.size())
                {
                    return false;
                }

                if (m_cycle_size == 0 && !prepare_cycle_cross())
                {
                    return false;
                }

                if (!m_in_second_region)
                {
                    if (m_j < m_a)
                    {
                        if (m_k < m_b)
                        {
                            m_current_element = m_cycle_offset + m_j * m_cycle_size + m_k;
                            ++m_k;
                            return true;
                        }
                        ++m_j;
                        m_k = m_a;
                        continue;
                    }
                    m_in_second_region = true;
                    m_j = m_a;
                    m_k = m_b;
                    continue;
                }

                if (m_j < m_b)
                {
                    if (m_k < m_cycle_size)
                    {
                        m_current_element = m_cycle_offset + m_j * m_cycle_size + m_k;
                        ++m_k;
                        return true;
                    }
                    ++m_j;
                    m_k = m_b;
                    continue;
                }

                ++m_cross_idx;
                m_cycle_size = 0;
            }
        }

        void advance()
        {
            if (m_is_end)
            {
                return;
            }

            if (m_phase == Phase::Cycle)
            {
                if (advance_cycle())
                {
                    return;
                }
                m_phase = Phase::End;
            }

            m_is_end = true;
        }

        std::reference_wrapper<const SetCoverCyc> m_owner;
        size_t m_set_index = 0;

        Phase m_phase = Phase::Cycle;
        bool m_is_end = false;
        size_t m_current_element = 0;

        size_t m_cross_idx = 0;
        bool m_in_second_region = false;
        size_t m_cycle_size = 0;
        size_t m_cycle_offset = 0;
        size_t m_a = 0;
        size_t m_b = 0;
        size_t m_j = 0;
        size_t m_k = 0;
    };

    SetCoverCyc(
        ull n_cycle_min_cuts,
        std::vector<std::vector<CycleCross<cycle_pos_T, cycle_id_T>>> cycle_crosses,
        std::vector<std::vector<cycle_pos_T>> cycle_positions,
        std::vector<size_t> cycle_sizes,
        std::vector<SetCost> set_weights)
        : m_nCycleMinCuts(n_cycle_min_cuts),
          m_cycleCrosses(std::move(cycle_crosses)),
          m_cyclePositions(std::move(cycle_positions)),
          m_cycleSizes(cycle_sizes),
          m_setWeights(std::move(set_weights)),
          m_cycleElementOffsets(build_cycle_element_offsets(this->m_cycleSizes)),
          m_cycleElementSpaceSize(compute_cycle_element_space_size(this->m_cycleSizes)) {
          };

private:
    const size_t m_nCycleMinCuts;
    const std::vector<std::vector<CycleCross<cycle_pos_T, cycle_id_T>>> m_cycleCrosses;
    const std::vector<std::vector<cycle_pos_T>> m_cyclePositions;
    const std::vector<size_t> m_cycleSizes;
    const std::vector<size_t> m_cycleElementOffsets;
    const size_t m_cycleElementSpaceSize;
    const std::vector<SetCost> m_setWeights;

public:
    const std::vector<size_t> &get_cycle_sizes() const
    {
        return m_cycleSizes;
    }
    const std::vector<std::vector<CycleCross<cycle_pos_T, cycle_id_T>>> &get_cycle_crosses() const
    {
        return m_cycleCrosses;
    }
    const std::vector<std::vector<cycle_pos_T>> &get_cycle_positions() const
    {
        return m_cyclePositions;
    }
    size_t get_num_sets() const noexcept
    {
        return m_setWeights.size();
    }
    size_t get_num_elements() const noexcept
    {
        return m_nCycleMinCuts;
    }
    SetCost get_set_cost(SetID set_index) const
    {
        return m_setWeights[set_index];
    }
    auto set_elements(size_t set_index) const
    {
        return std::ranges::subrange(SetElementsIterator(*this, set_index, false), SetElementsIterator(*this, set_index, true));
    }

    /// Apply a function to each covered cycle element in this set.
    /// Passes individual element indices for element-level operations on cycle coverage tracking.
    /// In this representation, there are more elements than there should be
    template <typename Callable>
    void for_each_covered_cycle_element(size_t set_index, Callable &&fn) const
    {
        // Iterate cycle crosses with (j, k) rectangles
        const auto &crosses = m_cycleCrosses[set_index];
        for (const auto &[cid, a, b] : crosses)
        {
            const size_t cycle_idx = static_cast<size_t>(cid);

            const size_t cycle_size = m_cycleSizes[cycle_idx];
            const size_t cycle_offset = m_cycleElementOffsets[cycle_idx];

            size_t a_clamped = static_cast<size_t>(a);
            size_t b_clamped = static_cast<size_t>(b);

            // Iterate first rectangle: [0, a_clamped) × [a_clamped, b_clamped)
            for (size_t j = 0; j < a_clamped; ++j)
            {
                for (size_t k = a_clamped; k < b_clamped; ++k)
                {
                    size_t element = cycle_offset + j * cycle_size + k;
                    fn(element);
                }
            }

            // Iterate second rectangle: [a_clamped, b_clamped) × [b_clamped, cycle_size)
            for (size_t j = a_clamped; j < b_clamped; ++j)
            {
                for (size_t k = b_clamped; k < cycle_size; ++k)
                {
                    size_t element = cycle_offset + j * cycle_size + k;
                    fn(element);
                }
            }
        }
    }

    void forEachElement(size_t set_index, const std::function<void(ElementID)> &func) const
    {
        for_each_covered_cycle_element(set_index, [&](size_t element)
                                       { func(element); });
    }

    size_t get_num_cycle_cuts() const
    {
        return m_nCycleMinCuts;
    }

private:
    static std::vector<size_t> build_cycle_element_offsets(const std::vector<size_t> &sizes)
    {
        std::vector<size_t> offsets(sizes.size(), 0);
        size_t prefix = 0;
        for (size_t i = 0; i < sizes.size(); ++i)
        {
            offsets[i] = prefix;
            prefix += sizes[i] * sizes[i];
        }
        return offsets;
    }

    static size_t compute_cycle_element_space_size(const std::vector<size_t> &sizes)
    {
        size_t total = 0;
        for (const auto sz : sizes)
        {
            total += sz * sz;
        }
        return total;
    }
};
