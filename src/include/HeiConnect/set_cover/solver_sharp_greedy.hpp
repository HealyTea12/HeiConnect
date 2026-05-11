#pragma once

#include "HeiConnect/set_cover/solver_base.hpp"
#include "HeiConnect/set_cover/common.hpp"

template <typename SetCoverType, typename Context, typename Solution>
concept GreedySharpContextCon = requires(Context context, size_t set_index, const SetCoverType &set_cover, Solution solution) {
    { context.is_element_covered(set_cover.get_num_elements() - 1) } -> std::convertible_to<bool>;
    { set_cover.get_num_elements() } -> std::convertible_to<size_t>;
    { set_cover.forEachSetCoveringElement(size_t{}, std::function<void(size_t)>{}) };
    { set_cover.get_set_cost(set_index) } -> std::convertible_to<double>;
    { solution.add_set(set_index) };
    { solution.get_solution_size() } -> std::convertible_to<size_t>;
};

// Context object could be passed to help it stop earlier for example
template <typename SetCoverType, typename Context, typename Solution>
class SetCoverSolverGreedySharp
{
public:
    struct Output
    {
        SolverStatus status = SolverStatus::Unknown;
    };
    Output solve(const SetCoverType &set_cover, const Context &context, Solution &solution)
    {
        Output output{};
        while (context.get_total_covered_elements() < set_cover.get_num_elements() &&
               solution.get_solution_size() != set_cover.get_num_sets())
        {
            for (size_t e = 0; e < set_cover.get_num_elements(); e++)
            {
                if (context.is_element_covered(e))
                    continue;
                auto cheapest_set_cost = std::numeric_limits<double>::max();
                size_t cheapest_set = std::numeric_limits<size_t>::max();
                set_cover.forEachSetCoveringElement(e, [&set_cover, &cheapest_set_cost, &cheapest_set](size_t set_index)
                                                    {
                    if (set_cover.get_set_cost(set_index) < cheapest_set_cost)
                    {
                        cheapest_set_cost = set_cover.get_set_cost(set_index);
                        cheapest_set = set_index;
                    } });
                if (cheapest_set == std::numeric_limits<size_t>::max())
                {
                    output.status = SolverStatus::Infeasible;
                    return output;
                }
                solution.add_set(cheapest_set);
                context.add_set(cheapest_set);
            }
        }
        output.status = SolverStatus::Feasible;
        return output;
    }
};
