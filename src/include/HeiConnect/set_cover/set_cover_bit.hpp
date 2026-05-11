#pragma once

#include "HeiConnect/set_cover/common.hpp"

#include <ranges>

struct SetCoverBit
{
    using ull = unsigned long long;
    using ElementID = size_t;
    using SetID = size_t;
    using SetCost = double;

    SetCoverBit(
        std::vector<ull> set_cover,
        size_t n_sets,
        size_t n_elements,
        std::vector<SetCost> costs)
        : m_adjacencyMatrix(std::move(set_cover)),
          m_nSets(n_sets),
          m_nElements(n_elements),
          m_nCols((n_elements - 1) / WORD_BITS + 1),
          m_costs(std::move(costs))
    {
    }

    size_t get_num_sets() const noexcept
    {
        return m_nSets;
    }

    size_t get_num_elements() const noexcept
    {
        return m_nElements;
    }

    SetCost get_set_cost(size_t set_index) const
    {
        return m_costs[set_index];
    }

    auto set_elements(size_t set_index)
    {
        return std::ranges::subrange(
            BitSetIterator(m_adjacencyMatrix.data() + set_index * m_nCols, 0, m_nCols),
            BitSetIterator(m_adjacencyMatrix.data() + set_index * m_nCols, 0, m_nCols, true));
    }

    BitSetIterator set_begin(size_t set_index)
    {
        return set_elements(set_index).begin();
    }

    BitSetIterator set_end(size_t set_index)
    {
        return set_elements(set_index).end();
    }

    void forEachElementBitMasked(size_t set_index, std::function<void(size_t)> func) const
    {
        for (size_t col = 0; col < m_nCols; col++)
        {
            func(m_adjacencyMatrix[set_index * m_nCols + col]);
        }
    }

    void forEachElement(size_t set_index, std::function<void(size_t)> func) const
    {
        for (size_t col = 0; col < m_nCols; col++)
        {
            ull word = m_adjacencyMatrix[set_index * m_nCols + col];
            const size_t base_element = col * WORD_BITS;
            while (word != 0)
            {
                size_t bit_idx = std::countr_zero(word);
                func(base_element + bit_idx);
                word ^= (1ULL << bit_idx);
            }
        }
    }

    ull get_col(size_t set_index, size_t col_index) const
    {
        return m_adjacencyMatrix[set_index * m_nCols + col_index];
    }

    size_t get_n_cols() const
    {
        return m_nCols;
    }

private:
    static constexpr size_t WORD_BITS = 8 * sizeof(ull);

    const std::vector<ull> m_adjacencyMatrix;
    const size_t m_nSets;
    const size_t m_nElements;
    const size_t m_nCols;
    const std::vector<SetCost> m_costs;
};
