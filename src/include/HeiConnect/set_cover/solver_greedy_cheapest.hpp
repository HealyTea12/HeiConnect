#pragma once

#include "HeiConnect/set_cover/solver_base.hpp"
#include "HeiConnect/set_cover/solver_greedy_context.hpp"
#include "HeiConnect/set_cover/trimmer.hpp"

#include <algorithm>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

template<typename SetCoverType, typename Context>
concept GreedyCheapestContextCon = requires(Context context, size_t set_index, const SetCoverType& set_cover) {
    { context.add_set(set_index) };
    { context.cover_count(set_index) } -> std::convertible_to<size_t>;
    { context.get_total_covered_elements() } -> std::convertible_to<size_t>;
    { set_cover.get_num_sets() } -> std::convertible_to<size_t>;
    { set_cover.get_set_cost(set_index) } -> std::convertible_to<double>;
    { set_cover.get_num_elements() } -> std::convertible_to<size_t>;
};

class SetCoverSolverGreedyCheapest
{
public:
    static constexpr std::string_view name = "Greedy cheapest solve";

    template<typename SetCoverType, typename Context>
        requires GreedyCheapestContextCon<SetCoverType, Context>
    auto operator()(std::shared_ptr<const SetCoverType> set_cover, Context context)
    {
        solve(*set_cover, context);
        return std::tuple{std::move(set_cover), std::move(context)};
    }

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
                const auto cost_a = set_cover.get_set_cost(a);
                const auto cost_b = set_cover.get_set_cost(b);
                return cost_a == cost_b ? a < b : cost_a < cost_b;
            });

        for (const auto set_index : set_indices)
        {
            if (context.get_total_covered_elements() == set_cover.get_num_elements())
            {
                break;
            }

            if (context.cover_count(set_index) > 0)
            {
                context.add_set(set_index);
            }
        }

        if (context.get_total_covered_elements() != set_cover.get_num_elements())
        {
            throw std::runtime_error("SetCoverSolverGreedyCheapest: set cover instance is infeasible.");
        }
    }
};
