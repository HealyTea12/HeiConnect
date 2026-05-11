#pragma once

#include "HeiConnect/set_cover/common.hpp"

#include <ranges>
#include <tuple>

template <class link_node_T, class link_edge_T>
    requires std::integral<link_node_T> && std::integral<link_edge_T>
class SetCoverPseudo
{
public:
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
        SetElementsIterator() = default;

        using value_type = size_t;
        using difference_type = std::ptrdiff_t;
        using iterator_concept = std::input_iterator_tag;
        using iterator_category = std::input_iterator_tag;

        SetElementsIterator(const SetCoverPseudo &owner, size_t set_index, bool is_end)
            : m_owner(owner), m_set_index(set_index), m_is_end(is_end)
        {
            if (!m_is_end)
            {
                load_word(0);
                next_bit();
            }
        }

        SetElementsIterator &operator++()
        {
            next_bit();
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
        static constexpr size_t WORD_BITS = 8 * sizeof(ull);

        void load_word(size_t word_index)
        {
            const auto &owner = m_owner.get();
            if (word_index >= owner.m_nCols)
            {
                m_is_end = true;
                return;
            }

            m_current_word_index = word_index;
            const auto &[u, v, w] = owner.m_links[m_set_index];
            (void)w;
            const size_t u_offset = static_cast<size_t>(u) * owner.m_nCols + m_current_word_index;
            const size_t v_offset = static_cast<size_t>(v) * owner.m_nCols + m_current_word_index;
            m_current_word = owner.m_partitionMatrix[u_offset] ^ owner.m_partitionMatrix[v_offset];

            const size_t remainder = owner.m_nMinCuts % WORD_BITS;
            if (m_current_word_index + 1 == owner.m_nCols && remainder != 0)
            {
                const ull mask = (1ULL << remainder) - 1ULL;
                m_current_word &= mask;
            }
        }

        void next_bit()
        {
            if (m_is_end)
            {
                return;
            }

            while (m_current_word == 0)
            {
                const size_t next_word = m_current_word_index + 1;
                load_word(next_word);
                if (m_is_end)
                {
                    return;
                }
            }

            const size_t bit = std::countr_zero(m_current_word);
            m_current_word ^= (1ULL << bit);
            m_current_element = m_current_word_index * WORD_BITS + bit;
        }

        std::reference_wrapper<const SetCoverPseudo> m_owner;
        size_t m_set_index = 0;
        size_t m_current_word_index = 0;
        size_t m_current_element = 0;
        ull m_current_word = 0;
        bool m_is_end = false;
    };

    static constexpr size_t WORD_BITS = 8 * sizeof(ull);
    SetCoverPseudo(
        const std::vector<ull> partiton_matrix,
        ull n_min_cuts,
        ull n_vertices,
        const std::vector<Link> links)
        : m_partitionMatrix(std::move(partiton_matrix)), m_nMinCuts(n_min_cuts),
          m_nVertices(n_vertices),
          m_links(std::move(links)), m_nCols(m_nMinCuts / WORD_BITS + (m_nMinCuts % WORD_BITS != 0)) {};
    size_t get_num_sets() const noexcept
    {
        return m_links.size();
    }
    size_t get_num_elements() const noexcept
    {
        return m_nMinCuts;
    }
    double get_set_cost(size_t set_index) const
    {
        return m_links[set_index].weight;
    }
    auto set_elements(size_t set_index) const
    {
        return std::ranges::subrange(SetElementsIterator(*this, set_index, false), SetElementsIterator(*this, set_index, true));
    }
    SetElementsIterator set_begin(size_t set_index) const
    {
        return SetElementsIterator(*this, set_index, false);
    }
    SetElementsIterator set_end(size_t set_index) const
    {
        return SetElementsIterator(*this, set_index, true);
    }
    void forEachElement(SetID set_index, const std::function<void(ElementID)> &func) const
    {
        auto [u, v, w] = m_links[set_index];
        (void)w;
        const size_t u_offset = static_cast<size_t>(u) * m_nCols;
        const size_t v_offset = static_cast<size_t>(v) * m_nCols;
        for (size_t word_index = 0; word_index < m_nCols; word_index++)
        {
            ull word = m_partitionMatrix[u_offset + word_index] ^ m_partitionMatrix[v_offset + word_index];
            size_t base_element = word_index * WORD_BITS;
            while (word != 0)
            {
                size_t bit_idx = std::countr_zero(word);
                func(base_element + bit_idx);
                word ^= (1ULL << bit_idx);
            }
        }
    }
    void forEachElementBitMasked(SetID set_index, const std::function<void(ElementID)> &func) const
    {
        auto [u, v, w] = m_links[set_index];
        (void)w;
        const size_t u_offset = static_cast<size_t>(u) * m_nCols;
        const size_t v_offset = static_cast<size_t>(v) * m_nCols;
        for (size_t word_index = 0; word_index < m_nCols; word_index++)
        {
            ull word = m_partitionMatrix[u_offset + word_index] ^ m_partitionMatrix[v_offset + word_index];
            func(word);
        }
    }
    ull get_col(size_t set_index, size_t col_index) const
    {
        auto [u, v, w] = m_links[set_index];
        (void)w;
        const size_t u_offset = static_cast<size_t>(u) * m_nCols + col_index;
        const size_t v_offset = static_cast<size_t>(v) * m_nCols + col_index;
        return m_partitionMatrix[u_offset] ^ m_partitionMatrix[v_offset];
    }
    size_t get_n_cols() const
    {
        return m_nCols;
    }
    const std::vector<Link> &get_links() const
    {
        return m_links;
    }
    const std::vector<ull> &get_partition_matrix() const
    {
        return m_partitionMatrix;
    }

private:
    const std::vector<ull> m_partitionMatrix;
    const ull m_nMinCuts;
    const size_t m_nCols;
    const size_t m_nVertices;
    const std::vector<Link> m_links;
};

static_assert(ForEachElementCon<SetCoverPseudo<uint64_t, uint64_t>>);
static_assert(ForEachElementBitMaskedCon<SetCoverPseudo<uint64_t, uint64_t>>);
