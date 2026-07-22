#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string_view>
#include <utility>

#include "HeiConnect/data_structures/immutable_graph.hpp"

namespace HeiConnect::visualization
{
    namespace detail
    {
        template<typename Graph>
        void write_edges(std::ostream& output, const Graph& graph, const std::string_view color)
        {
            std::set<std::pair<size_t, size_t>> edges;
            for (size_t u = 0; u < graph.num_vertices(); ++u)
            {
                for (size_t e = graph.graph.vertices[u]; e < graph.graph.vertices[u + 1]; ++e)
                {
                    const size_t v = graph.graph.edges[e];
                    edges.emplace(std::min(u, v), std::max(u, v));
                }
            }

            for (const auto& [u, v] : edges)
            {
                output << "    " << u << " -- " << v << " [color=" << color << "];\n";
            }
        }
    }

    template<typename OriginalGraph, typename LinkGraph>
    void write_to_dot(
        const OriginalGraph& original_graph,
        const LinkGraph& link_graph,
        const std::filesystem::path& output_path)
    {
        if (original_graph.num_vertices() != link_graph.num_vertices())
        {
            throw std::invalid_argument("The original graph and link graph must have the same number of vertices");
        }

        std::ofstream output{output_path};
        if (!output)
        {
            throw std::runtime_error("Could not open DOT output file: " + output_path.string());
        }

        output << "graph G {\n";
        for (size_t node = 0; node < original_graph.num_vertices(); ++node)
        {
            output << "    " << node << ";\n";
        }
        detail::write_edges(output, original_graph, "black");
        detail::write_edges(output, link_graph, "pink");
        output << "}\n";
    }
}
