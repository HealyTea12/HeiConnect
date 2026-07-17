#pragma once

#include <optional>
#include <string>
#include <tuple>

#include "HeiConnect/conn_aug/reducers/common.hpp"
#include "HeiConnect/conn_aug/reducers/full_mincut_dom_reducer.hpp"
#include "HeiConnect/conn_aug/reducers/full_reducer.hpp"
#include "HeiConnect/conn_aug/reducers/single_link_reducer.hpp"
#include "HeiConnect/data_structures/union_find.hpp"
#include "HeiConnect/pipeline/common.hpp"
#include "HeiConnect/set_cover/common.hpp"

template<int RecordStatsLevel = 0>
class FullSingleDomReducer
{
public:
    static std::string name()
    {
        return "Full Single Link Element Domination Reducer";
    }

    using LinkRemap = ConnAugLinkRemap;

    FullSingleDomReducer(bool project_in, bool project_out) : m_fullReducer(project_in, project_out)
    {}

    template<typename GraphType, typename LinkGraphType>
    auto operator()(const GraphType& graph, const LinkGraphType& link_graph)
    {
        LinkRemap link_remap = make_identity_link_remap(link_graph);
        UnionFind uf(graph.num_vertices());
        USSolution solution{};
        auto [reduced_graph, reduced_link_graph, ignored_link_remap, ignored_uf] =
            run(graph, link_graph, link_remap, uf, solution);
        (void)ignored_link_remap;
        (void)ignored_uf;
        return std::tuple{std::move(reduced_graph), std::move(reduced_link_graph)};
    }

    template<typename GraphType, typename LinkGraphType>
    auto operator()(
        const GraphType& graph,
        const LinkGraphType& link_graph,
        LinkRemap& link_remap,
        UnionFind& uf,
        USSolution& solution)
    {
        return run(graph, link_graph, link_remap, uf, solution);
    }

    template<typename NodeID, typename EdgeID, typename EdgeWeight, typename LinkEdgeID, typename LinkEdgeWeight>
    std::tuple<
        WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>,
        WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>,
        LinkRemap,
        UnionFind>
    run(const WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>& graph,
        const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph,
        LinkRemap& link_remap,
        UnionFind& uf,
        USSolution& solution)
    {
        if (link_remap.empty())
        {
            link_remap = make_identity_link_remap(link_graph);
        }

        if constexpr (RecordStatsLevel > 0)
        {
            m_totalRemovedByFull = 0;
            m_totalRemovedBySingle = 0;
            m_totalRemovedByDom = 0;
            m_metrics = StageMetrics{};
        }

        auto current_link_graph = link_graph;
        size_t iterations = 0;
        size_t num_changes = 0;

        while (true)
        {
            ++iterations;
            bool changed = false;

            const size_t before_full_edges = current_link_graph.num_edges();
            auto [full_graph, full_link_graph] = m_fullReducer.run(graph, current_link_graph);
            (void)full_graph;
            const size_t after_full_edges = full_link_graph.num_edges();
            if (!same_link_graph(current_link_graph, full_link_graph))
            {
                changed = true;
            }
            if constexpr (RecordStatsLevel > 0)
            {
                append_submetrics("full", iterations, m_fullReducer.emit_metrics());
            }
            current_link_graph = std::move(full_link_graph);

            const size_t before_single_edges = current_link_graph.num_edges();
            auto [single_graph, single_link_graph, single_link_remap, single_uf] =
                m_singleLinkReducer.run(graph, current_link_graph, link_remap, uf, solution);
            (void)single_graph;
            const size_t after_single_edges = single_link_graph.num_edges();
            link_remap = std::move(single_link_remap);
            uf = std::move(single_uf);
            if (!same_link_graph(current_link_graph, single_link_graph))
            {
                changed = true;
            }
            if constexpr (RecordStatsLevel > 0)
            {
                append_submetrics("single_link", iterations, m_singleLinkReducer.emit_metrics());
            }
            current_link_graph = std::move(single_link_graph);

            const size_t before_dom_edges = current_link_graph.num_edges();
            auto [dom_graph, dom_link_graph, dom_link_remap, dom_uf] =
                m_domReducer.run(graph, current_link_graph, link_remap, uf);
            (void)dom_graph;
            const size_t after_dom_edges = dom_link_graph.num_edges();
            link_remap = std::move(dom_link_remap);
            uf = std::move(dom_uf);
            if (!same_link_graph(current_link_graph, dom_link_graph))
            {
                changed = true;
            }
            if constexpr (RecordStatsLevel > 0)
            {
                append_submetrics("element_domination", iterations, m_domReducer.emit_metrics());
            }
            current_link_graph = std::move(dom_link_graph);

            if (changed)
            {
                ++num_changes;
            }

            if constexpr (RecordStatsLevel > 0)
            {
                m_totalRemovedByFull += before_full_edges - after_full_edges;
                m_totalRemovedBySingle += before_single_edges - after_single_edges;
                m_totalRemovedByDom += before_dom_edges - after_dom_edges;
            }

            if (!changed)
            {
                break;
            }
        }

        if constexpr (RecordStatsLevel > 0)
        {
            m_metrics->push_back({"num_iterations", std::to_string(iterations)});
            m_metrics->push_back({"num_changes", std::to_string(num_changes)});
            m_metrics->push_back(
                {"num_removed_links", std::to_string(link_graph.num_edges() - current_link_graph.num_edges())});
            m_metrics->push_back({"num_removed_by_full", std::to_string(m_totalRemovedByFull)});
            m_metrics->push_back({"num_removed_by_single_link", std::to_string(m_totalRemovedBySingle)});
            m_metrics->push_back({"num_removed_by_element_domination", std::to_string(m_totalRemovedByDom)});
        }

        return {graph, current_link_graph, link_remap, uf};
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
    void append_submetrics(const std::string& prefix, size_t iteration, const std::optional<StageMetrics>& metrics)
    {
        if (!metrics.has_value())
        {
            return;
        }
        for (const auto& metric : *metrics)
        {
            m_metrics->push_back(
                {prefix + "_iter_" + std::to_string(iteration) + "_" + metric.name, metric.printable_value});
        }
    }

    template<typename NodeID, typename LinkEdgeID, typename LinkEdgeWeight>
    bool same_link_graph(
        const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& lhs,
        const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& rhs) const
    {
        return lhs.graph.vertices == rhs.graph.vertices && lhs.graph.edges == rhs.graph.edges && lhs.weights == rhs.weights;
    }

    FullReducer<RecordStatsLevel> m_fullReducer;
    SingleLinkReducer<RecordStatsLevel> m_singleLinkReducer;
    FullMinCutDomReducer m_domReducer;
    std::optional<StageMetrics> m_metrics;
    size_t m_totalRemovedByFull{0};
    size_t m_totalRemovedBySingle{0};
    size_t m_totalRemovedByDom{0};
};
