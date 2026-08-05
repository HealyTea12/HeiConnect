#pragma once

#include "HeiConnect/set_cover/solver_base.hpp"
#include "HeiConnect/set_cover/solver_greedy_context.hpp"
#include "HeiConnect/set_cover/util.hpp"
#include "HeiConnect/pipeline/common.hpp"

#include <memory>
#include <queue>
#include <string_view>
#include <variant>

template<typename GreedyContext, typename SetCoverType>
concept GreedyContextCon = requires(GreedyContext context, size_t set_index) {
    { context.add_set(set_index) };
    { context.cover_count(set_index) } -> std::convertible_to<size_t>;
    { context.get_total_covered_elements() } -> std::convertible_to<size_t>;
};

template<size_t RecordMetricsLevel = 0>
class GreedySetCoverSolver
{
public:
    static constexpr std::string_view name = "Greedy solve";

    GreedySetCoverSolver() = default;

    template<typename SetCoverT, typename GreedyContext>
    auto operator()(std::shared_ptr<const SetCoverT> set_cover, GreedyContext context)
    {
        solve(*set_cover, context);
        return std::tuple{std::move(set_cover), std::move(context)};
    }

    std::optional<StageMetrics> emit_metrics() const
    {
        if constexpr (RecordMetricsLevel > 0)
        {
            return m_metrics;
        }
        else
        {
            return std::nullopt;
        }
    }

    template<typename SetCoverT, typename GreedyContext>
        requires GreedyContextCon<GreedyContext, SetCoverT>
    void solve(const SetCoverT& set_cover, GreedyContext& context)
    {
        using QueueEntry = std::pair<double, size_t>;
        std::vector<QueueEntry> initial_queue;
        initial_queue.reserve(set_cover.get_num_sets());
        for (size_t s{0}; s < set_cover.get_num_sets(); s++)
        {
            size_t covered = context.cover_count(s);
            const double cost_benefit_ratio = static_cast<double>(covered) / set_cover.get_set_cost(s);
            initial_queue.emplace_back(cost_benefit_ratio, s);
        }
        std::priority_queue<QueueEntry> pq{std::less<QueueEntry>{}, std::move(initial_queue)};

        while (context.get_total_covered_elements() < set_cover.get_num_elements() &&
               context.get_solution_size() != set_cover.get_num_sets())
        {
            auto best_set = pq.top();
            pq.pop();

            const auto best_set_idx = best_set.second;
            size_t covered = context.cover_count(best_set_idx);
            const double ratio = static_cast<double>(covered) / set_cover.get_set_cost(best_set_idx);
            if (ratio < best_set.first)
            {
                pq.push({ratio, best_set_idx});
                continue;
            }
            context.add_set(best_set_idx);
        }

        if constexpr (RecordMetricsLevel > 0)
        {
            m_metrics = StageMetrics{
                {"cost", std::to_string(HeiConnect::sc::cost(set_cover, context.get_solution()))},
                {"size", std::to_string(context.get_solution().size())}};
        }
    }

    template<typename SetCoverT, typename GreedyContext>
        requires GreedyContextCon<GreedyContext, SetCoverT>
    void repair(const SetCoverT& set_cover, GreedyContext& context)
    {
        solve(set_cover, context);
    }

private:
    std::optional<StageMetrics> m_metrics;
};
