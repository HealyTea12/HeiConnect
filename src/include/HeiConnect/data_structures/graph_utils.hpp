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
inline WeightedCRFGraph<> create_star_graph(size_t n_leaves)
{
    auto vertices = std::vector<size_t>(n_leaves + 2, 0);
    auto edges = std::vector<size_t>(2 * n_leaves, 0);
    auto weights = std::vector<double>(2 * n_leaves, 1.0);
    vertices[1] = n_leaves;
    for (size_t i = 1; i <= n_leaves; i++)
    {
        vertices[i + 1] = vertices[i] + 1;
        edges[i - 1] = i; // center to leaf
        edges[n_leaves + i - 1] = 0;
    }
    return WeightedCRFGraph<>{{vertices, edges}, weights};
}

inline WeightedCRFGraph<> create_random_tree(size_t n_nodes, unsigned int seed = 42)
{
    std::vector<bool> connected = std::vector<bool>(n_nodes, false);
    // tree has n - 1 edges
    auto vertices = std::vector<size_t>(n_nodes + 1, 0);
    auto edges = std::vector<size_t>(2 * (n_nodes - 1), 0);
    auto weights = std::vector<double>(2 * (n_nodes - 1), 1.0);
    auto random_engine = std::mt19937(seed);
    auto temp_edges = std::vector<std::vector<size_t>>(n_nodes, std::vector<size_t>{});
    std::uniform_int_distribution<size_t>
        dist{0, n_nodes - 1};
    for (size_t i = 0; i < n_nodes; i++)
    {
        // choose two nodes
        size_t u, v;
        do
        {
            u = dist(random_engine);
            v = dist(random_engine);
        } while (u == v || (connected[u] && connected[v]));
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