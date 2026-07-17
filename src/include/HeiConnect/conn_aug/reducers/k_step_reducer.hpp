#pragma once

#include "HeiConnect/data_structures/distance_oracle.hpp"
#include "HeiConnect/conn_aug/reducers/common.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"

/*
 THis does not look correct, but the idea is to take the k-step shortest distances
 between all pairs of nodes. If the k-step shortest distance d(u,v) is less than or equal to the link (u,v) weight,
 then we can remove the link (u,v)
*/
class KStepReducer
{
public:
    KStepReducer(size_t k) : m_k(k)
    {}
    template<typename GraphType, typename LinkGraphType>
    auto operator()(const GraphType& graph, const LinkGraphType& link_graph)
    {
        return run(graph, link_graph);
    }
    template<typename GraphType, typename LinkGraphType>
    auto run(const GraphType& graph, const LinkGraphType& link_graph)
    {
        using NodeID = typename GraphType::NodeID;
        using Distance = typename GraphType::EdgeWeight;
        // calculate all k-step shortest distances between all pairs of nodes in the link graph
        std::vector<Distance> k_step_distances = std::vector<Distance>(
            link_graph.num_vertices() * link_graph.num_vertices(),
            std::numeric_limits<Distance>::max());
        std::vector<bool> removable = std::vector<bool>(link_graph.num_edges(), false);
        for (NodeID u{0}; u < link_graph.num_vertices(); u++)
        {
            for (size_t step{0}; step < m_k; step++)
            {
                for (const auto& [v, w] : graph.get_neighbors(u))
                {
                    k_step_distances[u * link_graph.num_vertices() + v] =
                        std::min(k_step_distances[u * link_graph.num_vertices() + v], w);
                }
            }
        }
        for (NodeID u{0}; u < link_graph.num_vertices(); u++)
        {
            for (const auto& [v, w] : link_graph.get_neighbors(u))
            {
                if (k_step_distances[u * link_graph.num_vertices() + v] <= w)
                {
                    removable[link_graph.edge_id(u, v)] = true;
                }
            }
        }
        LinkGraphType new_link_graph = remove_links(link_graph, removable);
        return {graph, new_link_graph};
    }

private:
    size_t m_k;
};
