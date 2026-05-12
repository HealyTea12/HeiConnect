#pragma once

#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <vector>
#include <iostream>

template<bool Debug = false>
class ExhaustiveCombinationMoveGenerator
{
public:
    ExhaustiveCombinationMoveGenerator(
        std::size_t min_destroy_size = 1,
        std::size_t max_destroy_size = 0,
        std::size_t candidate_pool_size = 0,
        std::size_t max_moves = 0,
        bool expensive_first = true,
        bool include_full_solution_move = false) :
        m_minDestroySize(min_destroy_size),
        m_maxDestroySize(max_destroy_size),
        m_candidatePoolSize(candidate_pool_size),
        m_maxMoves(max_moves),
        m_expensiveFirst(expensive_first),
        m_includeFullSolutionMove(include_full_solution_move)
    {}

    template<typename SetCoverType, typename SolutionType, typename Callback>
    void forEachPotentialMove(const SetCoverType& set_cover, const SolutionType& solution, Callback&& callback) const
    {
        std::vector<std::size_t> selected_sets;
        selected_sets.reserve(solution.get_solution_size());

        for (auto set : solution.get_solution())
        {
            selected_sets.push_back(set);
        }

        sort_selected_sets(set_cover, selected_sets);

        if (m_candidatePoolSize > 0 && selected_sets.size() > m_candidatePoolSize)
        {
            selected_sets.resize(m_candidatePoolSize);
        }

        if (selected_sets.empty())
        {
            return;
        }

        const std::size_t selected_count = selected_sets.size();

        const std::size_t min_k = std::min(m_minDestroySize, selected_count);

        std::size_t max_k = (m_maxDestroySize == 0) ? selected_count : std::min(m_maxDestroySize, selected_count);

        if (!m_includeFullSolutionMove)
        {
            max_k = std::min(max_k, selected_count - 1);
        }

        if (min_k > max_k)
        {
            return;
        }

        std::size_t emitted_moves = 0;
        std::vector<std::size_t> current_move;
        current_move.reserve(max_k);

        for (std::size_t k = min_k; k <= max_k; ++k)
        {
            std::cout << "Emitting moves of size " << k << "..." << std::endl;
            const bool keep_going =
                enumerate_combinations_of_size(selected_sets, k, 0, current_move, emitted_moves, callback);

            if (!keep_going)
            {
                return;
            }

            if (move_limit_reached(emitted_moves))
            {
                return;
            }
        }
    }

private:
    template<typename SetCoverType>
    void sort_selected_sets(const SetCoverType& set_cover, std::vector<std::size_t>& selected_sets) const
    {
        if (m_expensiveFirst)
        {
            std::sort(selected_sets.begin(), selected_sets.end(), [&](std::size_t a, std::size_t b) {
                const auto cost_a = set_cover.get_set_cost(a);
                const auto cost_b = set_cover.get_set_cost(b);

                if (cost_a == cost_b)
                {
                    return a < b;
                }

                return cost_a > cost_b;
            });
        }
        else
        {
            std::sort(selected_sets.begin(), selected_sets.end());
        }
    }

    bool move_limit_reached(std::size_t emitted_moves) const
    {
        return m_maxMoves > 0 && emitted_moves >= m_maxMoves;
    }

    template<typename Callback>
    bool enumerate_combinations_of_size(
        const std::vector<size_t>& selected_sets,
        size_t target_size,
        size_t start_index,
        std::vector<size_t>& current_move,
        size_t& emitted_moves,
        Callback& callback) const
    {
        if (move_limit_reached(emitted_moves))
        {
            return false;
        }

        if (current_move.size() == target_size)
        {
            ++emitted_moves;

            // Convention:
            // callback(move) == true  -> continue enumeration
            // callback(move) == false -> stop enumeration
            return callback(current_move);
        }

        const size_t remaining_needed = target_size - current_move.size();

        const size_t last_possible_start = selected_sets.size() - remaining_needed;

        for (size_t i = start_index; i <= last_possible_start; ++i)
        {
            current_move.push_back(selected_sets[i]);

            const bool keep_going = enumerate_combinations_of_size(
                selected_sets,
                target_size,
                i + 1,
                current_move,
                emitted_moves,
                callback);

            current_move.pop_back();

            if (!keep_going)
            {
                return false;
            }
        }

        return true;
    }

    size_t m_minDestroySize;
    size_t m_maxDestroySize;
    size_t m_candidatePoolSize;
    size_t m_maxMoves;
    bool m_expensiveFirst;
    bool m_includeFullSolutionMove;
};