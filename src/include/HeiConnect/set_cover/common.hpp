#pragma once

#include <algorithm>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <immintrin.h>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <omp.h>

#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/minmax.hpp"
#include "HeiConnect/sc_reduction/cactus_min_cuts.hpp"

struct SolverConfig
{
    std::size_t seed = 0;
    std::size_t time_limit_seconds = 0;
    std::size_t thread_count = 1;
};

template<typename T>
concept SetCoverCon = requires(T t, size_t i) {
    { t.get_num_sets() } -> std::convertible_to<size_t>;
    { t.get_num_elements() } -> std::convertible_to<size_t>;
    { t.get_set_cost(i) } -> std::floating_point;
};

template<typename ContextType>
concept IncrementallyRemovableContext = requires(ContextType& context, size_t set_index) {
    context.remove_set(set_index);
};

template<typename ContextType, typename SolutionType>
concept CanRemoveContext = requires(
    const ContextType& context,
    size_t set_index,
    const SolutionType& solution) {
    { context.can_remove(set_index, solution) } -> std::convertible_to<bool>;
};

class USSolution
{
public:
    using SetID = size_t;

    void add_set(SetID set_index)
    {
        m_solution.insert(set_index);
    }

    void remove_set(SetID set_index)
    {
        m_solution.erase(set_index);
    }

    const std::unordered_set<SetID>& get_solution() const noexcept
    {
        return m_solution;
    }

    size_t get_solution_size() const noexcept
    {
        return m_solution.size();
    }

    size_t size() const noexcept
    {
        return m_solution.size();
    }

    void reset()
    {
        m_solution.clear();
    }

private:
    std::unordered_set<SetID> m_solution;
};


class VectorSolution
{
public:
    using SetID = size_t;

    void add_set(SetID set_index)
    {
        m_solution.emplace_back(set_index);
    }

    void remove_set(SetID set_index)
    {
        m_solution.erase(std::remove(m_solution.begin(), m_solution.end(), set_index), m_solution.end());
    }

    const std::vector<SetID>& get_solution() const noexcept
    {
        return m_solution;
    }

    size_t get_solution_size() const noexcept
    {
        return m_solution.size();
    }

    size_t size() const noexcept
    {
        return m_solution.size();
    }

    void reset()
    {
        m_solution.clear();
    }

private:
    std::vector<SetID> m_solution;
};


template<typename ContextType, typename SolutionType = USSolution>
class BoundContext
{
public:
    BoundContext(ContextType context, SolutionType solution) :
        m_context(std::move(context)),
        m_solution(std::move(solution))
    {}

    ContextType& get_context() noexcept
    {
        return m_context;
    }

    const ContextType& get_context() const noexcept
    {
        return m_context;
    }


    // Forward context operations
    void add_set(size_t set_index)
    {
        m_context.add_set(set_index);
        m_solution.add_set(set_index);
    }

    void remove_set(size_t set_index)
        requires IncrementallyRemovableContext<ContextType>
    {
        m_context.remove_set(set_index);
        m_solution.remove_set(set_index);
    }

    template<typename SolType>
        requires CanRemoveContext<ContextType, SolType>
    bool can_remove(size_t set_index, const SolType& solution) const
    {
        return m_context.can_remove(set_index, solution);
    }

    bool can_remove(size_t set_index) const
        requires CanRemoveContext<ContextType, SolutionType>
    {
        return m_context.can_remove(set_index, m_solution);
    }

    // Remove a selected set without changing a context that cannot remove sets incrementally.
    void remove_set_from_solution(size_t set_index)
    {
        m_solution.remove_set(set_index);
    }

    void rebuild_context()
        requires requires(ContextType& context, size_t set_index) {
            context.reset();
            context.add_set(set_index);
        }
    {
        m_context.reset();
        for (const size_t set_index : m_solution.get_solution())
        {
            m_context.add_set(set_index);
        }
    }


    using SolutionContainer = decltype(std::declval<SolutionType>().get_solution());

    SolutionContainer& get_solution()
    {
        return m_solution.get_solution();
    }

    const SolutionContainer& get_solution() const
    {
        return m_solution.get_solution();
    }

    size_t get_solution_size() const
    {
        return m_solution.get_solution_size();
    }

    BoundContext& operator=(const BoundContext& other)
    {
        if (this != &other)
        {
            m_solution = other.m_solution;
        }
        return *this;
    }

    size_t cover_count(size_t set_index) const
    {
        return m_context.cover_count(set_index);
    }

    size_t get_total_covered_elements() const
    {
        return m_context.get_total_covered_elements();
    }

    bool is_element_covered(size_t element_index) const
    {
        return m_context.is_element_covered(element_index);
    }

    void reset()
    {
        m_context.reset();
        m_solution.reset();
    }

private:
    ContextType m_context;
    SolutionType m_solution;
};

template<typename SetCoverType>
concept ForEachElementCon = requires(SetCoverType t, size_t s) {
    { t.forEachElement(s, std::function<void(size_t)>{}) };
};

template<typename SetCoverType>
concept ForEachElementBitMaskedCon = requires(SetCoverType t, size_t s) {
    { t.forEachElementBitMasked(s, std::function<void(size_t)>{}) };
};

template<typename SolverT, typename SetCoverType>
concept SetCoverSolverCon = requires(SolverT solver, const SetCoverType& sc, const SolverConfig& config) {
    { solver.solve(sc, config) };
};

enum class SolverStatus
{
    Optimal,
    Feasible,
    Infeasible,
    TimeLimit,
    Error,
    Unknown
};


template<typename Range, typename T>
concept range_of = std::ranges::input_range<Range> && std::same_as<std::ranges::range_value_t<Range>, T>;

using ull = unsigned long long;

template<typename NewRangeElementType, typename Range>
auto discretize_weights(const Range& weights, size_t num_bins)
{
    auto max_weight = *std::ranges::max_element(weights);
    return weights | std::views::transform([num_bins, max_weight](auto w) {
               return static_cast<NewRangeElementType>(std::ceil(w / max_weight * num_bins));
           });
}

class BitSetIterator
{
public:
    using value_type = size_t;
    using difference_type = std::ptrdiff_t;
    using iterator_concept = std::input_iterator_tag;
    using iterator_category = std::input_iterator_tag;

    BitSetIterator() = default;

    BitSetIterator(const ull* data, size_t word_index, size_t n_words, bool is_end = false) :
        m_data(data),
        current_word_idx(word_index),
        current_word(0),
        m_is_end(is_end),
        m_n_words(n_words)
    {
        if (current_word_idx < m_n_words)
        {
            current_word = m_data[current_word_idx];
        }
        else
        {
            m_is_end = true;
        }
        next_bit();
    };
    BitSetIterator& operator++()
    {
        next_bit();
        return *this;
    }
    BitSetIterator operator++(int)
    {
        BitSetIterator tmp = *this;
        ++(*this);
        return tmp;
    }
    size_t operator*() const
    {
        return current_word_idx * sizeof(ull) * 8 + current_bit_idx;
    }
    bool operator!=(const BitSetIterator& other) const
    {
        return m_is_end != other.m_is_end;
    }
    friend bool operator==(const BitSetIterator& a, const BitSetIterator& b)
    {
        return !(a != b);
    }

private:
    void next_bit()
    {
        if (m_is_end)
            return;
        while (current_word == 0 && !m_is_end)
        {
            current_word_idx++;
            if (current_word_idx >= m_n_words)
            {
                m_is_end = true;
                return;
            }
            current_word = m_data[current_word_idx];
        }
        current_bit_idx = std::countr_zero(current_word);
        current_word ^= (1ULL << current_bit_idx);
    }
    size_t current_word_idx;
    ull current_word;
    size_t current_bit_idx = 0;
    const ull* m_data;
    bool m_is_end;
    size_t m_n_words;
};

static_assert(std::input_iterator<BitSetIterator>);

template<typename SetCoverType>
concept ElementIterableCon = requires(SetCoverType t, size_t s) {
    { t.element_begin(s) } -> std::input_iterator;
    { t.element_end(s) } -> std::input_iterator;
    { t.set_elements(s) } -> std::ranges::input_range;
};
