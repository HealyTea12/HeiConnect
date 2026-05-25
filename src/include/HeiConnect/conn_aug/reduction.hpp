
#pragma once

#include <algorithm>
#include <vector>

#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/data_structures/distance_oracle.hpp"

template<int RecordStatsLevel = 0>
class BasicLinkDomReducer
{
public:
    template<typename NodeID, typename EdgeID, typename EdgeWeight, typename LinkEdgeID, typename LinkEdgeWeight>
    std::tuple<WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>, WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>>
    run(const WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>& graph,
        const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph)
    {
        TableDistOracle distance_oracle(graph.graph);
        std::vector<std::tuple<NodeID, NodeID, LinkEdgeWeight>> link_edges = link_graph.csr_to_vec_links();
        std::sort(link_edges.begin(), link_edges.end(), [&](const auto& a, const auto& b) {
            auto [u_a, v_a, w_a] = a;
            auto [u_b, v_b, w_b] = b;
            return w_a < w_b;
        });
        std::vector<bool> removable = std::vector<bool>(link_edges.size(), false);
        for (size_t i{0}; i < link_edges.size(); i++)
        {
            auto [u_i, v_i, w_i] = link_edges[i];
            for (size_t j{0}; j < i; j++)
            {
                auto [u_j, v_j, w_j] = link_edges[j];
                if (dominates(link_edges[j], link_edges[i], distance_oracle))
                {
                    removable[i] = true;
                    break;
                }
            }
        }
        if constexpr (RecordStatsLevel > 0)
        {
            size_t num_removed = std::count(removable.begin(), removable.end(), true);
            m_stats.num_removed_links = num_removed;
        }
        auto new_link_graph = remove_links<NodeID, LinkEdgeID, LinkEdgeWeight>(link_graph, removable);
        return {graph, new_link_graph};
    }


private:
    template<typename NodeID, typename LinkEdgeWeight>
    bool dominates(
        const std::tuple<NodeID, NodeID, LinkEdgeWeight>& link_i,
        const std::tuple<NodeID, NodeID, LinkEdgeWeight>& link_j,
        const TableDistOracle& distance_oracle)
    {
        auto [u_i, v_i, w_i] = link_i;
        auto [u_j, v_j, w_j] = link_j;
        return on_path(u_i, v_i, u_j, distance_oracle) && on_path(u_i, v_i, v_j, distance_oracle);
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
    struct Stats
    {
        size_t num_removed_links{0};
    } m_stats;
};
