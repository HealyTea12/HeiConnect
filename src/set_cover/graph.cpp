#include <unordered_set>
#include <random>
#include "set_cover/graph.hpp"

/// Generate a graph consisting of all non-existing edges in the input graph.
/// The weights of the new edges are determined by the provided edge_weight function.
/// @param graph
/// @param edge_weight
/// @return
WeightedCRFGraph generate_links(const CRFGraph &graph, double (*edge_weight)(size_t u, size_t v))
{
    std::vector<size_t> new_vertices{};
    std::vector<size_t> new_edges{};
    std::vector<double> weights{};
    new_vertices.push_back(0);
    for (size_t u = 0; u < graph.vertices.size() - 1; ++u)
    {
        std::unordered_set<size_t> neighbours{};
        for (size_t idx = graph.vertices[u]; idx < graph.vertices[u + 1]; ++idx)
        {
            neighbours.insert(graph.edges[idx]);
        }
        for (size_t i = 0; i < graph.vertices.size() - 1; i++)
        {
            if (i == u)
                continue;
            if (neighbours.find(i) == neighbours.end())
            {
                auto weight = edge_weight(u, i);
                weights.push_back(weight);
                new_edges.push_back(i);
            }
        }
        new_vertices.push_back(new_edges.size());
    }
    return WeightedCRFGraph{{new_vertices, new_edges}, weights};
}

WeightedCRFGraph generate_links(const WeightedCRFGraph &graph, double (*edge_weight)(size_t u, size_t v))
{
    return generate_links(graph.graph, edge_weight);
}

WeightedCRFGraph add_links(const WeightedCRFGraph &graph,
                           const WeightedCRFGraph &link_graph,
                           const std::unordered_set<size_t> &selected_edges)
{
    auto new_vertices = std::vector<size_t>(graph.graph.vertices.size(), 0);
    std::vector<size_t> new_edges{};
    std::vector<double> new_weights{};
    for (size_t u = 0; u < graph.graph.vertices.size() - 1; ++u)
    {
        for (size_t e = graph.graph.vertices[u]; e < graph.graph.vertices[u + 1]; ++e)
        {
            auto v = graph.graph.edges[e];
            new_edges.emplace_back(v);
            new_weights.emplace_back(graph.weights[e]);
        }
        for (size_t e = link_graph.graph.vertices[u]; e < link_graph.graph.vertices[u + 1]; ++e)
        {
            if (selected_edges.find(e) != selected_edges.end())
            {
                auto v = link_graph.graph.edges[e];
                new_edges.emplace_back(v);
                new_weights.emplace_back(0.0); // link weights are just the cost of adding them
            }
        }
        new_vertices[u + 1] = new_edges.size();
    }
    return WeightedCRFGraph{
        {new_vertices, new_edges},
        new_weights};
}

// Make the graph bidirectional by adding reverse edges
WeightedCRFGraph make_bidirectional(const WeightedCRFGraph &graph)
{
    auto new_vertices = std::vector<size_t>(graph.graph.vertices.size(), 0);
    std::vector<size_t> new_edges{};
    std::vector<double> new_weights{};
    for (size_t u = 0; u < graph.graph.vertices.size() - 1; ++u)
    {
        // add original edges
        for (size_t e = graph.graph.vertices[u]; e < graph.graph.vertices[u + 1]; ++e)
        {
            auto v = graph.graph.edges[e];
            new_edges.emplace_back(v);
            new_weights.emplace_back(graph.weights[e]);
        }
        // add reverse edges
        for (size_t v = 0; v < graph.graph.vertices.size() - 1; ++v)
        {
            for (size_t e = graph.graph.vertices[v]; e < graph.graph.vertices[v + 1]; ++e)
            {
                if (graph.graph.edges[e] == u)
                {
                    new_edges.emplace_back(v);
                    new_weights.emplace_back(graph.weights[e]);
                }
            }
        }
        new_vertices[u + 1] = new_edges.size();
    }
    return WeightedCRFGraph{
        {new_vertices, new_edges},
        new_weights};
}