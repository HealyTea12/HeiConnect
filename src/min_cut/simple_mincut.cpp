#include "simple_mincut.hpp"
#include <limits>

double global_mincut_simple(const WeightedCRFGraph &graph)
{
    // nodes that have not yet been contracted
    std::vector<bool> active_nodes(graph.graph.vertices.size() - 1, true);
    auto total_active_nodes = graph.graph.vertices.size() - 1;
    auto min_cut = std::numeric_limits<double>::max();
    // store added weight to edges during contraction
    std::vector<double> extra_weights(graph.graph.edges.size(), 0.0);
    while (total_active_nodes > 1)
    {
        // Phase of the Stoer-Wagner algorithm
        std::vector<double> weights(graph.graph.vertices.size() - 1, 0.0);
        // set A as named in the paper
        std::vector<bool> in_a_set(graph.graph.vertices.size() - 1, false);
        size_t prev = 0;
        size_t last = 0;
        size_t next = 0;
        for (auto k = 0ull; k < graph.graph.vertices.size() - 1; ++k)
        {
            if (active_nodes[k])
            {
                next = k;
                in_a_set[k] = true;
                for (size_t e = graph.graph.vertices[k]; e < graph.graph.vertices[k + 1]; ++e)
                {
                    size_t v = graph.graph.edges[e];
                    if (active_nodes[v] && !in_a_set[v]) // maybe don't need this check
                    {
                        weights[v] += graph.weights[e] + extra_weights[e];
                    }
                }
                break;
            }
        }
        for (size_t i = 0; i < total_active_nodes - 1; ++i)
        {
            double max_weight = -1.0;
            // set the starting node to the first node still in the graph
            for (size_t j = 0; j < graph.graph.vertices.size() - 1; ++j)
            {
                if (active_nodes[j] && !in_a_set[j] && weights[j] > max_weight)
                {
                    max_weight = weights[j];
                    next = j;
                }
            }

            in_a_set[next] = true;
            if (i == total_active_nodes - 2) // last iteration
            {
                // Update min cut
                if (weights[next] < min_cut)
                {
                    min_cut = weights[next];
                }

                // Merge prev and next
                for (size_t e = graph.graph.vertices[next]; e < graph.graph.vertices[next + 1]; ++e)
                {
                    size_t v = graph.graph.edges[e];
                    for (size_t f = graph.graph.vertices[prev]; f < graph.graph.vertices[prev + 1]; ++f)
                    {
                        if (graph.graph.edges[f] == v)
                        {
                            extra_weights[f] += graph.weights[e];
                            for (size_t h = graph.graph.vertices[v]; h < graph.graph.vertices[v + 1]; ++h)
                            {
                                if (graph.graph.edges[h] == prev)
                                {
                                    extra_weights[h] += graph.weights[e];
                                    break;
                                }
                            }
                            break;
                        }
                    }
                }
                active_nodes[next] = false;
                total_active_nodes--;
            }
            else
            {
                for (size_t e = graph.graph.vertices[next]; e < graph.graph.vertices[next + 1]; ++e)
                {
                    size_t v = graph.graph.edges[e];
                    if (active_nodes[v] && !in_a_set[v])
                    {
                        weights[v] += graph.weights[e] + extra_weights[e];
                    }
                }
                prev = next;
            }
        }
    }

    return min_cut;
}