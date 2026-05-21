#pragma once

namespace HeiConnect::sc
{
    template<typename SetCoverType, typename SolutionType>
    typename SetCoverType::SetCost cost(const SetCoverType& set_cover, const SolutionType& solution)
    {
        auto cost = typename SetCoverType::SetCost{};
        if constexpr (requires { solution.get_solution(); })
        {
            for (size_t set : solution.get_solution())
            {
                cost += set_cover.get_set_cost(set);
            }
        }
        else
        {
            for (size_t set : solution)
            {
                cost += set_cover.get_set_cost(set);
            }
        }
        return cost;
    }
}