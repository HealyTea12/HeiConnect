#pragma once

#include "HeiConnect/set_cover/common.hpp"
#include "HeiConnect/set_cover/trimmer_context.hpp"
// Default trimmer context uses the same basic context as the greedy solver.
#include "HeiConnect/set_cover/solver_greedy_context.hpp"

#include <algorithm>
#include <unordered_set>
#include <vector>

template <typename SetCoverT, typename Solution, typename TrimmerContext>
concept TrimmerRequirements = requires(const SetCoverT &set_cover, Solution &solution, TrimmerContext &context) {
    { set_cover.get_num_sets() } -> std::convertible_to<size_t>;
    { set_cover.get_num_elements() } -> std::convertible_to<size_t>;
    { set_cover.get_set_cost(size_t{}) } -> std::floating_point;
    { set_cover.forEachElement(size_t{}, std::function<void(size_t)>{}) };
    { solution.get_solution() } -> std::ranges::input_range;
    { solution.remove_set(size_t{}) };
    { context.add_set(size_t{}) };
    { context.remove_set(size_t{}) };
    { context.can_remove(size_t{}, solution) } -> std::convertible_to<bool>;
};

template <typename SetCoverT, typename Solution, typename TrimmerContext = BasicContext<SetCoverT, Solution>>
    requires TrimmerRequirements<SetCoverT, Solution, TrimmerContext>
class SetCoverTrimmer
{
public:
    static void trim(const SetCoverT &set_cover, Solution &solution)
    {
        TrimmerContext context{std::make_shared<const SetCoverT>(set_cover)};
        trim(set_cover, solution, context);
    }

    static void populate_and_trim(const SetCoverT &set_cover, Solution &solution, TrimmerContext &context)
    {
        std::vector<size_t> selected_sets(solution.get_solution().begin(), solution.get_solution().end());

        for (size_t set_index : selected_sets)
        {
            context.add_set(set_index);
        }
        trim(set_cover, solution, context);
    }

    static void trim(const SetCoverT &set_cover, Solution &solution, TrimmerContext &context)
    {
        std::vector<size_t> selected_sets(solution.get_solution().begin(), solution.get_solution().end());

        std::sort(selected_sets.begin(), selected_sets.end(), [&](size_t a, size_t b)
                  { return set_cover.get_set_cost(a) > set_cover.get_set_cost(b); });

        for (size_t set_index : selected_sets)
        {
            if (context.can_remove(set_index, solution))
            {
                context.remove_set(set_index);
                solution.remove_set(set_index);
            }
        }
    }
};
