#pragma once

#include "HeiConnect/data_structures/immutable_graph.hpp"


template<typename NodeID, typename LinkEdgeID, typename LinkEdgeWeight>
WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>
remove_links(const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph, const std::vector<bool>& removable)
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