
#pragma once

#include <algorithm>
#include <vector>

#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/data_structures/distance_oracle.hpp"

template<int RecordStatsLevel = 0>
class BasicLinkDomReducer
{
    template<typename NodeID, typename EdgeID, typename EdgeWeight, typename LinkEdgeID, typename LinkEdgeWeight>
    void
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
                if (dominates(link_edges[i], link_edges[j], distance_oracle))
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

public:
    struct Stats
    {
        size_t num_removed_links{0};
    } m_stats;
};
