#pragma once

#include <chrono>
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
#include "HeiConnect/tools/timer.hpp"

/**
 * TODO: FOrced links are contracted one at a time. There is definitely a more efficient way of doing it i.e:
 * not reconstructing the whole graph each time, just building the union find for all of them and then
 * materializing the contractions all at once.
 * The single-link phase now builds one union find and materializes all selected links in bulk.
 */
template<int RecordStatsLevel = 0>
class FullSingleDomReducer
{
public:
    static std::string name()
    {
        return "Full Single Link Element Domination Reducer";
    }

    using LinkRemap = ConnAugLinkRemap;

    FullSingleDomReducer() = default;

    explicit FullSingleDomReducer(ConnectivityAugmentationReductionConfig config) :
        m_config(config),
        m_fullReducer(config)
    {}

    FullSingleDomReducer(bool project_in, bool project_out) : m_fullReducer(project_in, project_out)
    {
        m_config.run_project_in = project_in;
        m_config.run_project_out = project_out;
    }

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
        using Clock = std::chrono::high_resolution_clock;
        auto setup_start = Clock::now();
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

        auto current_graph = graph;
        auto current_link_graph = link_graph;
        uf = UnionFind(current_graph.num_vertices());
        auto setup_end = Clock::now();
        size_t iterations = 0;
        size_t num_changes = 0;

        if constexpr (RecordStatsLevel > 0)
        {
            m_metrics->push_back(
                {"setup_time",
                 HeiConnect::tools::format_duration(setup_start, setup_end, HeiConnect::tools::TimeUnit::Seconds)});
        }

        while (true)
        {
            if (m_config.max_rounds > 0 && iterations >= m_config.max_rounds)
            {
                break;
            }

            ++iterations;
            bool changed = false;

            const size_t before_full_edges = current_link_graph.num_edges();
            auto full_start = Clock::now();
            auto [full_graph, full_link_graph] = m_fullReducer.run(current_graph, current_link_graph);
            auto full_end = Clock::now();
            (void)full_graph;
            const size_t after_full_edges = full_link_graph.num_edges();
            if (!same_link_graph(current_link_graph, full_link_graph))
            {
                changed = true;
            }
            if constexpr (RecordStatsLevel > 0)
            {
                append_submetrics("full", iterations, m_fullReducer.emit_metrics());
                append_time_metric("full", iterations, full_start, full_end);
            }
            current_link_graph = std::move(full_link_graph);

            const size_t before_single_edges = current_link_graph.num_edges();
            auto single_start = Clock::now();
            auto single_reducer_start = Clock::now();
            std::vector<LinkEdgeID> selected_links;
            if (m_config.run_single_link)
            {
                selected_links = m_singleLinkReducer.run(current_graph, current_link_graph);
            }
            auto single_reducer_end = Clock::now();
            auto single_materialize_start = Clock::now();
            if (!selected_links.empty())
            {
                const auto links = current_link_graph.csr_to_vec_links();
                for (const LinkEdgeID link_id : selected_links)
                {
                    const auto& [u, v, weight] = links[link_id];
                    (void)weight;
                    solution.add_set(link_remap.at(normalize_link(u, v)).original_id);
                }

                UnionFind selection_uf(current_graph.num_vertices());
                add_links_to_union_find_frozen_stack(
                    current_graph,
                    current_link_graph,
                    selection_uf,
                    selected_links);
                auto [single_contracted_graph, single_contracted_link_graph, single_node_remap] =
                    materialize_contractions(current_graph, current_link_graph, selection_uf, link_remap);
                (void)single_node_remap;
                changed = true;
                current_graph = std::move(single_contracted_graph);
                current_link_graph = std::move(single_contracted_link_graph);
                uf = UnionFind(current_graph.num_vertices());
            }
            auto single_materialize_end = Clock::now();
            const size_t after_single_edges = current_link_graph.num_edges();
            auto single_end = Clock::now();
            if constexpr (RecordStatsLevel > 0)
            {
                if (m_config.run_single_link)
                {
                    append_submetrics("single_link", iterations, m_singleLinkReducer.emit_metrics());
                }
                append_time_metric(
                    "single_link_reducer",
                    iterations,
                    single_reducer_start,
                    single_reducer_end);
                append_time_metric(
                    "single_link_materialize",
                    iterations,
                    single_materialize_start,
                    single_materialize_end);
                append_time_metric("single_link", iterations, single_start, single_end);
            }

            const size_t before_dom_edges = current_link_graph.num_edges();
            auto element_domination_start = Clock::now();
            if (m_config.run_element_domination)
            {
                auto [dom_graph, dom_link_graph, dom_link_remap, dom_uf] =
                    m_domReducer.run(current_graph, current_link_graph, link_remap, uf);
                (void)dom_graph;
                link_remap = std::move(dom_link_remap);
                uf = std::move(dom_uf);
                auto [contracted_graph, contracted_link_graph, node_remap] =
                    materialize_contractions(current_graph, dom_link_graph, uf, link_remap);
                (void)node_remap;
                if (!same_link_graph(current_link_graph, contracted_link_graph))
                {
                    changed = true;
                }
                current_graph = std::move(contracted_graph);
                current_link_graph = std::move(contracted_link_graph);
                uf = UnionFind(current_graph.num_vertices());
            }
            auto element_domination_end = Clock::now();
            const size_t after_dom_edges = current_link_graph.num_edges();
            if constexpr (RecordStatsLevel > 0)
            {
                if (m_config.run_element_domination)
                {
                    append_submetrics("element_domination", iterations, m_domReducer.emit_metrics());
                }
                append_time_metric("element_domination", iterations, element_domination_start, element_domination_end);
            }

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

        return {current_graph, current_link_graph, link_remap, uf};
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
    template<typename TimePoint>
    void append_time_metric(const std::string& prefix, size_t iteration, TimePoint start, TimePoint end)
    {
        m_metrics->push_back(
            {prefix + "_iter_" + std::to_string(iteration) + "_total_time",
             HeiConnect::tools::format_duration(start, end, HeiConnect::tools::TimeUnit::Seconds)});
    }

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
        return lhs.graph.vertices == rhs.graph.vertices && lhs.graph.edges == rhs.graph.edges &&
            lhs.weights == rhs.weights;
    }

    ConnectivityAugmentationReductionConfig m_config{};
    FullReducer<RecordStatsLevel> m_fullReducer;
    SingleLinkReducer<RecordStatsLevel> m_singleLinkReducer;
    FullMinCutDomReducer m_domReducer;
    std::optional<StageMetrics> m_metrics;
    size_t m_totalRemovedByFull{0};
    size_t m_totalRemovedBySingle{0};
    size_t m_totalRemovedByDom{0};
};
