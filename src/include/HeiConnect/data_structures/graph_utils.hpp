#pragma once
#include "data_structures/immutable_graph.hpp"

WeightedCRFGraph<> create_cycle_graph(size_t n_nodes)
{
    auto vertices = std::vector<size_t>(n_nodes, 0);
    std::vector<size_t> edges{};
    std::vector<double> weights{};
    for (size_t u = 0; u < n_nodes; u++)
    {
        edges.push_back((u + 1) % n_nodes);
        weights.push_back(1.0);
        vertices[u + 1] = edges.size();
    }
    return WeightedCRFGraph<>{CRFGraph<>{vertices, edges}, weights};
}
