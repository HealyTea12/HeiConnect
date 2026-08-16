
#pragma once

#include <chrono>
#include <cstdint>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include "HeiConnect/conn_aug/reducers/common.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/pipeline/common.hpp"
#include "HeiConnect/set_cover/common.hpp"
#include "HeiConnect/tools/timer.hpp"
#include "HeiConnect/data_structures/union_find.hpp"
#include "HeiConnect/sc_reduction/transform_single_builders.hpp"


/*  This reducer verifies whether a minimum cut/element is covered by a single link/set.
    If so, it adds that link/set to the solution.

 *  INITIAL BASIC IMPLEMNTATION
 */
template<int RecordStatsLevel = 0>
class SingleLinkReducer
{
public:
    static std::string name()
    {
        return "Single Link Reducer";
    }

    SingleLinkReducer()
    {}

    template<typename GraphType, typename LinkGraphType>
    auto operator()(const GraphType& graph, const LinkGraphType& link_graph)
    {
        return run(graph, link_graph);
    }

    template<typename NodeID, typename EdgeID, typename EdgeWeight, typename LinkEdgeID, typename LinkEdgeWeight>
    std::vector<LinkEdgeID>
    run(const WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>& graph,
        const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph)
    {
        using Clock = std::chrono::high_resolution_clock;

        uint64_t n_selected_links{};
        auto start_set_cover_oracle = Clock::now();
        auto sc = construct_set_cover_oracle(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        auto end_set_cover_oracle = Clock::now();

        auto start_select_links = Clock::now();
        std::vector<size_t> covering_sets(sc.get_num_sets());
        std::vector<bool> is_selected_link(link_graph.num_edges(), false);
        std::vector<LinkEdgeID> selected_links{};
        for (size_t element_id{}; element_id < sc.get_num_elements(); element_id++)
        {
            covering_sets.clear();
            sc.forEachSet(element_id, [&](size_t set_id) { covering_sets.emplace_back(set_id); });
            if (covering_sets.size() == 1)
            {
                is_selected_link[covering_sets[0]] = true;
            }
        }
        for (size_t link_id{}; link_id < is_selected_link.size(); ++link_id)
        {
            if (is_selected_link[link_id])
            {
                ++n_selected_links;
                selected_links.emplace_back(static_cast<LinkEdgeID>(link_id));
            }
        }
        auto end_select_links = Clock::now();

        if constexpr (RecordStatsLevel > 0)
        {
            m_metrics = StageMetrics{
                {"num_selected_links", std::to_string(n_selected_links)},
                {
                    "time_set_cover_oracle",
                    HeiConnect::tools::format_duration(
                        start_set_cover_oracle,
                        end_set_cover_oracle,
                        HeiConnect::tools::TimeUnit::Seconds),
                },
                {
                    "time_select_links",
                    HeiConnect::tools::format_duration(
                        start_select_links,
                        end_select_links,
                        HeiConnect::tools::TimeUnit::Seconds),
                },
            };
        }

        return selected_links;
    }

    std::optional<StageMetrics> emit_metrics() const
    {
        if constexpr (RecordStatsLevel > 0)
        {
            return m_metrics;
        }
        else
        {
            return std::nullopt;
        }
    }

private:
    std::optional<StageMetrics> m_metrics;
};
