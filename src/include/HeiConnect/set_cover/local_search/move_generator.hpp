#pragma once

#include <algorithm>
#include <vector>
#include <cstddef>

#include "small_set.hpp"

class SingleSetRemovalMoveGenerator
{
public:
    explicit SingleSetRemovalMoveGenerator(bool expensive_first = true) : m_expensiveFirst(expensive_first)
    {}

    template<typename SetCoverType, typename SolutionType, typename Callback>
    void forEachPotentialMove(const SetCoverType& set_cover, const SolutionType& solution, Callback&& callback) const
    {
        std::vector<std::size_t> sets;
        sets.reserve(solution.get_solution_size());

        for (auto set : solution.get_solution())
        {
            sets.push_back(set);
        }

        if (m_expensiveFirst)
        {
            std::sort(sets.begin(), sets.end(), [&](std::size_t a, std::size_t b) {
                return set_cover.get_set_cost(a) > set_cover.get_set_cost(b);
            });
        }
        else
        {
            std::sort(sets.begin(), sets.end());
        }

        for (auto set : sets)
        {
            SmallMove<1> move;
            move.push_back(set);

            if (!callback(move))
            {
                return;
            }
        }
    }

private:
    bool m_expensiveFirst;
};