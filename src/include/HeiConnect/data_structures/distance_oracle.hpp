#pragma once

#include "HeiConnect/data_structures/immutable_graph.hpp"

class TableDistOracle
{
    using Distance = uint64_t;

public:
    TableDistOracle(const CRFGraph<size_t, size_t>& graph)
    {
        std::vector<std::tuple<size_t, size_t>> edges = graph.edges_vector();
        n = graph.vertices.size();
        m_distances.resize(n * n, std::numeric_limits<Distance>::max());
        for (size_t i = 0; i < n; ++i)
        {
            m_distances[i * n + i] = 0;
        }
        for (const auto& edge : edges)
        {
            auto [u, v] = edge;
            m_distances[u * n + v] = 1;
            m_distances[v * n + u] = 1;
        }
        // Floyd-Warshall algorithm to compute all-pairs shortest paths
        for (size_t k = 0; k < n; ++k)
        {
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < n; ++j)
                {
                    if (m_distances[i * n + j] > m_distances[i * n + k] + m_distances[k * n + j])
                    {
                        m_distances[i * n + j] = m_distances[i * n + k] + m_distances[k * n + j];
                    }
                }
            }
        }
    }
    Distance get_distance(size_t u, size_t v) const
    {
        return m_distances[u * n + v];
    };

private:
    std::vector<Distance> m_distances;
    int n;
};
