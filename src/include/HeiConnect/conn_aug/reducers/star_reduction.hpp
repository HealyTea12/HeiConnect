#pragma once
#include <algorithm>
#include <tuple>
#include <vector>
#include <queue>
#include <array>
#include <limits>
#include <string_view>
#include <optional>

#include "HeiConnect/data_structures/immutable_graph.hpp"

#include "common.hpp"

template<int RecordStatsLevel = 0>
class StarReducer
{
public:
    static constexpr std::string_view name = "Star reduction";

    StarReducer(size_t degree_threshold = 10) : m_degreeThreshold(degree_threshold)
    {}

    template<typename GraphType, typename LinkGraphType>
    auto operator()(const GraphType& graph, const LinkGraphType& link_graph) const
    {
        return run(graph, link_graph);
    }

    template<typename NodeID, typename EdgeID, typename EdgeWeight, typename LinkEdgeID, typename LinkEdgeWeight>
    std::tuple<WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>, WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>>
    run(const WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>& graph,
        const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph) const
    {
        std::vector<bool> removable_edges = std::vector<bool>(link_graph.num_edges(), false);
        // identify high degree vertices
        std::vector<NodeID> high_degree_vertices{};
        high_degree_vertices.reserve(graph.num_vertices());
        size_t highest_degree = 0;
        NodeID highest_degree_vertex = 0;
        for (NodeID u{0}; u < graph.num_vertices(); ++u)
        {
            if (graph.degree(u) > m_degreeThreshold)
            {
                high_degree_vertices.emplace_back(u);
            }
            if (graph.degree(u) > highest_degree)
            {
                highest_degree = graph.degree(u);
                highest_degree_vertex = u;
            }
        }
        // Do BFS from highest degree node and remove all edges to other high degree node
        std::vector<EdgeWeight> min_weight =
            std::vector<EdgeWeight>(graph.num_vertices(), std::numeric_limits<EdgeWeight>::max());
        for (NodeID u : high_degree_vertices)
        {
            auto [parent, depth] = bfs_parent_depth(graph.graph, u);
            calculate_min_weight_star(graph, link_graph, u, parent, depth, min_weight);
            for (NodeID i{0}; i < link_graph.num_vertices(); i++)
            {
                for (LinkEdgeID e{link_graph.graph.vertices[i]}; e < link_graph.graph.vertices[i + 1]; e++)
                {
                    NodeID j = link_graph.graph.edges[e];
                    if (depth[i] == 1 && depth[j] == 1)
                    {
                        if (min_weight[i] + min_weight[j] < link_graph.weights[e])
                        {
                            removable_edges[e] = true;
                        }
                    }
                }
            }
        }
        auto new_link_graph = remove_links<NodeID, LinkEdgeID, LinkEdgeWeight>(link_graph, removable_edges);

        if constexpr (RecordStatsLevel > 0)
        {
            size_t num_removed = std::count(removable_edges.begin(), removable_edges.end(), true);
            m_metrics = StageMetrics{{"num_removed_links", std::to_string(num_removed)}};
        }
        return {graph, new_link_graph};
    }

private:
    template<typename NodeID>
    std::pair<NodeID, NodeID>
    projected_edge(NodeID i, NodeID j, const std::vector<NodeID>& parent, const std::vector<size_t>& depth, NodeID u)
        const
    {
        std::array<NodeID, 2> endpoints;
        int k = 0;
        while (i != j)
        {
            if (depth[i] >= depth[j])
            {
                if (i == u)
                {
                    endpoints[k++] = parent[i];
                }
                else if (parent[i] == u)
                {
                    endpoints[k++] = i;
                }
                i = parent[i];
            }
            else
            {
                if (j == u)
                {
                    endpoints[k++] = parent[j];
                }
                else if (parent[j] == u)
                {
                    endpoints[k++] = j;
                }
                j = parent[j];
            }
        }
        if (k == 2)
        {
            return {endpoints[0], endpoints[1]};
        }
        else
        {
            return {std::numeric_limits<NodeID>::max(), std::numeric_limits<NodeID>::max()};
        }
    };
    template<typename NodeID, typename EdgeID, typename EdgeWeight, typename LinkEdgeID, typename LinkEdgeWeight>
    void calculate_min_weight_star(
        const WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>& graph,
        const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph,
        NodeID u,
        const std::vector<NodeID>& parent,
        const std::vector<size_t>& depth,
        std::vector<EdgeWeight>& min_weight) const
    {
        for (NodeID i{0}; i < link_graph.num_vertices(); i++)
        {
            for (LinkEdgeID e{link_graph.graph.vertices[i]}; e < link_graph.graph.vertices[i + 1]; e++)
            {
                NodeID j = link_graph.graph.edges[e];
                LinkEdgeWeight w = link_graph.weights[e];
                auto [a, b] = projected_edge(i, j, parent, depth, u);
                if (a != std::numeric_limits<NodeID>::max())
                {
                    min_weight[a] = std::min(min_weight[a], w);
                    min_weight[b] = std::min(min_weight[b], w);
                }
            }
        }
    }

    template<typename NodeID, typename EdgeID>
    std::tuple<std::vector<NodeID>, std::vector<size_t>>
    bfs_parent_depth(const CRFGraph<NodeID, EdgeID>& graph, NodeID start) const
    {
        std::vector<NodeID> parent(graph.vertices.size() - 1, std::numeric_limits<NodeID>::max());
        std::vector<size_t> depth(graph.vertices.size() - 1, std::numeric_limits<size_t>::max());
        std::queue<NodeID> q;
        parent[start] = start;
        depth[start] = 0;
        q.push(start);
        while (!q.empty())
        {
            NodeID u = q.front();
            q.pop();
            for (EdgeID e{graph.vertices[u]}; e < graph.vertices[u + 1]; ++e)
            {
                NodeID v = graph.edges[e];
                if (parent[v] == std::numeric_limits<NodeID>::max())
                {
                    parent[v] = u;
                    depth[v] = depth[u] + 1;
                    q.push(v);
                }
            }
        }
        return {parent, depth};
    }

public:
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
    size_t m_degreeThreshold;
    mutable std::optional<StageMetrics> m_metrics;
};