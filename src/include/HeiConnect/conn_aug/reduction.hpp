
#pragma once

#include <algorithm>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/data_structures/distance_oracle.hpp"

template<int RecordStatsLevel = 0>
class BasicLinkDomReducer
{
public:
    static constexpr std::string_view name = "Data pruning";

    template<typename GraphType, typename LinkGraphType>
    auto operator()(const GraphType& graph, const LinkGraphType& link_graph)
    {
        return run(graph, link_graph);
    }

    template<typename NodeID, typename EdgeID, typename EdgeWeight, typename LinkEdgeID, typename LinkEdgeWeight>
    std::tuple<WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>, WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>>
    run(const WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>& graph,
        const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph)
    {
        TableDistOracle distance_oracle(graph.graph);
        struct LinkRecord
        {
            size_t original_index;
            NodeID u;
            NodeID v;
            LinkEdgeWeight w;
        };

        std::vector<LinkRecord> link_edges;
        link_edges.reserve(link_graph.num_edges());
        auto original_link_edges = link_graph.csr_to_vec_links();
        for (size_t i = 0; i < original_link_edges.size(); ++i)
        {
            auto [u, v, w] = original_link_edges[i];
            link_edges.push_back({i, u, v, w});
        }
        std::sort(link_edges.begin(), link_edges.end(), [&](const auto& a, const auto& b) { return a.w < b.w; });
        std::vector<bool> removable = std::vector<bool>(link_graph.num_edges(), false);
        for (size_t i{0}; i < link_edges.size(); i++)
        {
            const auto& link_i = link_edges[i];
            for (size_t j{0}; j < i; j++)
            {
                const auto& link_j = link_edges[j];
                if (dominates(link_j, link_i, distance_oracle))
                {
                    removable[link_i.original_index] = true;
                    break;
                }
            }
        }
        if constexpr (RecordStatsLevel > 0)
        {
            size_t num_removed = std::count(removable.begin(), removable.end(), true);
            m_metrics = StageMetrics{{"num_removed_links", std::to_string(num_removed)}};
        }
        auto new_link_graph = remove_links<NodeID, LinkEdgeID, LinkEdgeWeight>(link_graph, removable);
        return {graph, new_link_graph};
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
    template<typename LinkRecord>
    bool dominates(const LinkRecord& link_i, const LinkRecord& link_j, const TableDistOracle& distance_oracle)
    {
        return on_path(link_i.u, link_i.v, link_j.u, distance_oracle) &&
            on_path(link_i.u, link_i.v, link_j.v, distance_oracle);
    }

    template<typename NodeID>
    bool on_path(NodeID u, NodeID v, NodeID x, const TableDistOracle& distance_oracle)
    {
        return distance_oracle.get_distance(u, x) + distance_oracle.get_distance(x, v) ==
            distance_oracle.get_distance(u, v);
    }


    template<typename NodeID, typename LinkEdgeID, typename LinkEdgeWeight>
    WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight> remove_links(
        const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph,
        const std::vector<bool>& removable)
    {
        std::vector<NodeID> new_vertices(link_graph.num_vertices() + 1, 0);
        std::vector<LinkEdgeID> new_edges{};
        new_edges.reserve(link_graph.num_edges());
        std::vector<LinkEdgeWeight> new_weights{};
        new_edges.reserve(link_graph.num_edges());
        for (NodeID u{0}; u < link_graph.num_vertices(); u++)
        {
            for (LinkEdgeID e{link_graph.graph.vertices[u]}; e < link_graph.graph.vertices[u + 1]; e++)
            {
                if (!removable[e])
                {
                    new_edges.push_back(link_graph.graph.edges[e]);
                    new_weights.push_back(link_graph.weights[e]);
                }
            }
            new_vertices[u + 1] = new_edges.size();
        }
        return {{new_vertices, new_edges}, new_weights};
    }

public:
    std::optional<StageMetrics> m_metrics;
};
