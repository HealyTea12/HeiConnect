#pragma once

#include "HeiConnect/pipeline/common.hpp"
#include "HeiConnect/set_cover/common.hpp"
#include "HeiConnect/set_cover/set_cover_cyc.hpp"
#include "HeiConnect/set_cover/set_cover_double.hpp"
#include "HeiConnect/set_cover/set_cover_pseudo.hpp"
#include "HeiConnect/set_cover/trimmer_context.hpp"
#include "HeiConnect/set_cover/solver_greedy_context.hpp"
#include "HeiConnect/set_cover/util.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <variant>
#include <unordered_set>
#include <vector>

template<typename SetCoverT, typename TrimmerContext>
concept TrimmerRequirements = requires(const SetCoverT& set_cover, TrimmerContext& context) {
    { set_cover.get_num_sets() } -> std::convertible_to<size_t>;
    { set_cover.get_num_elements() } -> std::convertible_to<size_t>;
    { set_cover.get_set_cost(size_t{}) } -> std::floating_point;
    { set_cover.forEachElement(size_t{}, std::function<void(size_t)>{}) };
    { context.add_set(size_t{}) };
    { context.remove_set(size_t{}) };
    { context.can_remove(size_t{}) } -> std::convertible_to<bool>;
};

template<size_t RecordMetricsLevel = 0>
class SetCoverTrimmer
{
    template<typename LinkNode, typename LinkEdge, typename LinkWeight>
    using DoubleSetCover = SetCoverDouble<
        SetCoverPseudo<LinkNode, LinkEdge>,
        SetCoverCyc<LinkNode, LinkEdge, LinkWeight>>;

    template<typename LinkNode, typename LinkEdge, typename LinkWeight>
    using DoubleContext = ContextDouble<
        SetCoverPseudo<LinkNode, LinkEdge>,
        SetCoverCyc<LinkNode, LinkEdge, LinkWeight>,
        BitPackedContext<SetCoverPseudo<LinkNode, LinkEdge>>,
        CycContext<LinkNode, LinkEdge, LinkWeight>>;

    template<typename LinkNode, typename LinkEdge, typename LinkWeight>
    using BoundDoubleContext = BoundContext<DoubleContext<LinkNode, LinkEdge, LinkWeight>, USSolution>;

public:
    static constexpr std::string_view name = "Trim redundant sets";
    SetCoverTrimmer() = default;

    template<typename SetCoverT, typename TrimmerContext>
    auto operator()(std::shared_ptr<const SetCoverT> set_cover, TrimmerContext context)
    {
        trim(*set_cover, context);
        return std::tuple{std::move(set_cover), std::move(context)};
    }

    template<typename SetCoverT, typename TrimmerContext>
        requires TrimmerRequirements<SetCoverT, TrimmerContext>
    void trim(const SetCoverT& set_cover, TrimmerContext& context)
    {
        std::vector<size_t> selected_sets(context.get_solution().begin(), context.get_solution().end());

        std::sort(selected_sets.begin(), selected_sets.end(), [&](size_t a, size_t b) {
            const auto cost_a = set_cover.get_set_cost(a);
            const auto cost_b = set_cover.get_set_cost(b);
            return cost_a == cost_b ? a < b : cost_a > cost_b;
        });

        for (size_t set_index : selected_sets)
        {
            if (context.can_remove(set_index))
            {
                context.remove_set(set_index);
            }
        }

        record_metrics(set_cover, context);
    }

    template<typename LinkNode, typename LinkEdge, typename LinkWeight>
    void trim(
        const DoubleSetCover<LinkNode, LinkEdge, LinkWeight>& set_cover,
        BoundDoubleContext<LinkNode, LinkEdge, LinkWeight>& context)
    {
        std::vector<size_t> selected_sets(context.get_solution().begin(), context.get_solution().end());

        std::sort(selected_sets.begin(), selected_sets.end(), [&](size_t a, size_t b) {
            const auto cost_a = set_cover.get_set_cost(a);
            const auto cost_b = set_cover.get_set_cost(b);
            return cost_a == cost_b ? a < b : cost_a > cost_b;
        });

        for (const size_t set_index : selected_sets)
        {
            if (covers_without(set_cover, context.get_solution(), set_index))
            {
                context.remove_set(set_index);
            }
        }

        auto rebuilt_context = make_double_context(set_cover);
        for (const size_t set_index : context.get_solution())
        {
            rebuilt_context.add_set(set_index);
        }
        context.get_context() = std::move(rebuilt_context);

        record_metrics(set_cover, context);
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

private:
    template<typename LinkNode, typename LinkEdge, typename LinkWeight>
    static DoubleContext<LinkNode, LinkEdge, LinkWeight>
    make_double_context(const DoubleSetCover<LinkNode, LinkEdge, LinkWeight>& set_cover)
    {
        BitPackedContext pseudo_context{set_cover.first()};
        CycContext cycle_context{set_cover.second()};
        return DoubleContext<LinkNode, LinkEdge, LinkWeight>{
            set_cover,
            std::move(pseudo_context),
            std::move(cycle_context)};
    }

    template<typename LinkNode, typename LinkEdge, typename LinkWeight>
    static bool covers_without(
        const DoubleSetCover<LinkNode, LinkEdge, LinkWeight>& set_cover,
        const std::unordered_set<size_t>& solution,
        size_t omitted_set)
    {
        auto trial_context = make_double_context(set_cover);
        for (const size_t set_index : solution)
        {
            if (set_index == omitted_set)
            {
                continue;
            }

            trial_context.add_set(set_index);
            if (trial_context.get_total_covered_elements() == set_cover.get_num_elements())
            {
                return true;
            }
        }
        return trial_context.get_total_covered_elements() == set_cover.get_num_elements();
    }

    template<typename SetCoverT, typename TrimmerContext>
    void record_metrics(const SetCoverT& set_cover, const TrimmerContext& context)
    {
        if constexpr (RecordMetricsLevel > 0)
        {
            m_metrics = StageMetrics{
                {"cost", std::to_string(HeiConnect::sc::cost(set_cover, context.get_solution()))},
                {"size", std::to_string(context.get_solution().size())}
            };
        }
    }
    std::optional<StageMetrics> m_metrics;
};

template<size_t RecordMetricsLevel = 0>
SetCoverTrimmer() -> SetCoverTrimmer<RecordMetricsLevel>;
