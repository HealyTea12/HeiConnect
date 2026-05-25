#pragma once

#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/bfs.hpp"

class TableDistOracle
{
public:
    using Distance = uint64_t;

    TableDistOracle(const CRFGraph<size_t, size_t>& graph)
    {
        std::vector<std::tuple<size_t, size_t>> edges = graph.edges_vector();
        n = graph.vertices.size() - 1;
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
        // BFS from each vertex: set distances when a node is discovered
        std::vector<bool> visited{};
        visited.resize(n, false);
        std::vector<size_t> queue{};
        queue.reserve(n);

        for (size_t i = 0; i < n; ++i)
        {
            std::fill(visited.begin(), visited.end(), false);
            visited[i] = true;
            queue.clear();
            queue.push_back(i);
            while (!queue.empty())
            {
                size_t u = queue.front();
                queue.erase(queue.begin());
                for (size_t e = graph.vertices[u]; e < graph.vertices[u + 1]; ++e)
                {
                    size_t v = graph.edges[e];
                    if (!visited[v])
                    {
                        visited[v] = true;
                        m_distances[i * n + v] = m_distances[i * n + u] + 1;
                        m_distances[v * n + i] = m_distances[i * n + v];
                        queue.push_back(v);
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
    size_t n{0};
};
