#pragma once

#include "HeiConnect/set_cover/solver_base.hpp"
#include "HeiConnect/set_cover/solver_greedy_context.hpp"
#include "HeiConnect/set_cover/trimmer.hpp"

#include <memory>

template<typename SetCoverType, typename Context>
concept GreedyCheapestContextCon = requires(Context context, size_t set_index, const SetCoverType& set_cover) {
    { context.add_set(set_index) };
    { context.get_total_covered_elements() } -> std::convertible_to<size_t>;
    { set_cover.get_num_sets() } -> std::convertible_to<size_t>;
    { set_cover.get_set_cost(set_index) } -> std::convertible_to<double>;
    { set_cover.get_num_elements() } -> std::convertible_to<size_t>;
    { context.get_solution_size() } -> std::convertible_to<size_t>;
};

class SetCoverSolverGreedyCheapest
{
public:
    template<typename SetCoverType, typename Context>
        requires GreedyCheapestContextCon<SetCoverType, Context>
    void solve(const SetCoverType& set_cover, Context& context)
    {
        std::vector<typename SetCoverType::SetID> set_indices =
            std::vector<typename SetCoverType::SetID>(set_cover.get_num_sets(), 0);
        std::iota(set_indices.begin(), set_indices.end(), 0);
        std::sort(
            set_indices.begin(),
            set_indices.end(),
            [&](typename SetCoverType::SetID a, typename SetCoverType::SetID b) {
                return set_cover.get_set_cost(a) < set_cover.get_set_cost(b);
            });
        size_t i = 0;
        while (context.get_total_covered_elements() < set_cover.get_num_elements() &&
               context.get_solution_size() < set_cover.get_num_sets())
        {
            context.add_set(set_indices[i]);
            i++;
        }
    }
};
