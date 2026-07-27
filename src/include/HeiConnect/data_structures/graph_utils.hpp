#pragma once
#include <random>

#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/graph.hpp"

inline WeightedCRFGraph<> create_cycle_graph(size_t n_nodes)
{
    auto vertices = std::vector<size_t>(n_nodes + 1, 0);
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

inline WeightedCRFGraph<> create_cycle_graph_undirected(size_t n_nodes)
{
    auto vertices = std::vector<size_t>(n_nodes + 1, 0);
    std::vector<size_t> edges{};
    std::vector<double> weights{};
    for (size_t u = 0; u < n_nodes; u++)
    {
        edges.push_back((u + 1) % n_nodes);
        weights.push_back(1.0);
        edges.push_back((u + n_nodes - 1) % n_nodes);
        weights.push_back(1.0);
        vertices[u + 1] = edges.size();
    }
    return WeightedCRFGraph<>{CRFGraph<>{vertices, edges}, weights};
}

// Creates a star graph with one center node (node 0) connected to n_leaves leaf nodes.
template <class node_T = uint64_t, class edge_T = uint64_t, class weight_T = double>
inline WeightedCRFGraph<node_T, edge_T, weight_T> create_star_graph(node_T n_leaves)
{
    auto vertices = std::vector<edge_T>(n_leaves + 2, static_cast<edge_T>(0));
    auto edges = std::vector<node_T>(2 * n_leaves, static_cast<node_T>(0));
    auto weights = std::vector<weight_T>(2 * n_leaves, static_cast<weight_T>(1.0));
    vertices[1] = n_leaves;
    for (node_T i = static_cast<node_T>(1); i <= n_leaves; i++)
    {
        vertices[i + 1] = vertices[i] + 1;
        edges[i - 1] = i; // center to leaf
        edges[n_leaves + i - 1] = 0;
    }
    return WeightedCRFGraph<node_T, edge_T, weight_T>{{vertices, edges}, weights};
}

inline WeightedCRFGraph<> create_random_tree(size_t n_nodes, unsigned int seed = 42)
{
    // tree has n - 1 edges
    auto vertices = std::vector<size_t>(n_nodes + 1, 0);
    auto edges = std::vector<size_t>(2 * (n_nodes - 1), 0);
    auto weights = std::vector<double>(2 * (n_nodes - 1), 1.0);
    auto random_engine = std::mt19937(seed);
    auto temp_edges = std::vector<std::vector<size_t>>(n_nodes, std::vector<size_t>{});
    for (size_t u = 1; u < n_nodes; u++)
    {
        // choose one of the already connected nodes
        std::uniform_int_distribution<size_t> dist{0, u - 1};
        auto v = dist(random_engine);
        temp_edges[u].emplace_back(v);
        temp_edges[v].emplace_back(u);
    }
    for (size_t u = 0; u < n_nodes; u++)
    {
        auto starting_point = vertices[u]; // we start placing us neighbours here
        vertices[u + 1] = starting_point + temp_edges[u].size();
        for (size_t i = 0; i < temp_edges[u].size(); i++)
        {
            auto v = temp_edges[u][i];
            edges[starting_point + i] = v;
        }
    }
    return WeightedCRFGraph<>{{vertices, edges}, weights};
}
