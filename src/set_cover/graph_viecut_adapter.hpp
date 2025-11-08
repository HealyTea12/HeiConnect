#pragma once
#include "extern/VieCut/lib/data_structure/mutable_graph.h"
#include "set_cover/graph.hpp"

mutable_graph create_mutable_graph(const CRFGraph &crf_graph)
{
    mutable_graph g{};
    g.start_construction(crf_graph.vertices.size() - 1, crf_graph.edges.size());
    for (size_t u = 0; u < crf_graph.vertices.size() - 1; ++u)
    {
        g.new_node();
        for (size_t e = crf_graph.vertices[u]; e < crf_graph.vertices[u + 1]; ++e)
        {
            size_t v = crf_graph.edges[e];
            g.new_edge(u, v, 1.0);
        }
    }
    g.finish_construction();
    return g;
}

mutable_graph create_mutable_graph(const WeightedCRFGraph &graph)
{
    mutable_graph g{};
    g.start_construction(graph.graph.vertices.size() - 1, graph.graph.edges.size());
    for (size_t u = 0; u < graph.graph.vertices.size() - 1; ++u)
    {
        g.new_node();
        for (size_t e = graph.graph.vertices[u]; e < graph.graph.vertices[u + 1]; ++e)
        {
            size_t v = graph.graph.edges[e];
            g.new_edge(u, v, graph.weights[e]);
        }
    }
    g.finish_construction();
    return g;
}