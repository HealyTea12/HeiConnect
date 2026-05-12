#pragma once
/*
#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "small_set.hpp"

template<std::size_t MaxK>
class TopKDestroyMoveGenerator
{
public:
    TopKDestroyMoveGenerator(std::size_t max_destroy_size, std::size_t candidate_pool_size, std::size_t max_moves) :
        m_maxDestroySize(max_destroy_size),
        m_candidatePoolSize(candidate_pool_size),
        m_maxMoves(max_moves)
    {}

    template<typename SetCoverType, typename SolutionType, typename Callback>
    void forEachPotentialMove(const SetCoverType& set_cover, const SolutionType& solution, Callback&& callback) const
    {
        std::vector<std::size_t> selected;
        selected.reserve(solution.get_solution_size());

        for (auto set : solution.get_solution())
        {
            selected.push_back(set);
        }

        std::sort(selected.begin(), selected.end(), [&](std::size_t a, std::size_t b) {
            return set_cover.get_set_cost(a) > set_cover.get_set_cost(b);
        });

        if (selected.size() > m_candidatePoolSize)
        {
            selected.resize(m_candidatePoolSize);
        }

        std::size_t emitted = 0;

        DestroyMove<MaxK> move;

        for (std::size_t k = 1; k <= m_maxDestroySize && k <= MaxK; ++k)
        {
            if (!enumerate_combinations(selected, 0, k, move, emitted, callback))
            {
                return;
            }

            if (emitted >= m_maxMoves)
            {
                return;
            }
        }
    }

private:
    template<typename Callback>
    bool enumerate_combinations(
        const std::vector<std::size_t>& selected,
        std::size_t start,
        std::size_t target_size,
        DestroyMove<MaxK>& move,
        std::size_t& emitted,
        Callback&& callback) const
    {
        if (emitted >= m_maxMoves)
        {
            return false;
        }

        if (move.size == target_size)
        {
            ++emitted;
            return callback(move);
        }

        for (std::size_t i = start; i < selected.size(); ++i)
        {
            move.removed[move.size++] = selected[i];

            const bool keep_going = enumerate_combinations(selected, i + 1, target_size, move, emitted, callback);

            --move.size;

            if (!keep_going)
            {
                return false;
            }
        }

        return true;
    }

    std::size_t m_maxDestroySize;
    std::size_t m_candidatePoolSize;
    std::size_t m_maxMoves;
};
*/