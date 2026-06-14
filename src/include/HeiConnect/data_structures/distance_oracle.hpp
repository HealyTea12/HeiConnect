#pragma once

#include "HeiConnect/data_structures/immutable_graph.hpp"

#include <limits>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>


/*
 * A simple distance oracle that stores all pairwise distances in a table.
 * The constructor computes the distances using BFS from each vertex.
 */
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

template<typename NodeID = size_t, typename EdgeID = size_t, typename WeightType = double>
class WeightedTableDistOracle
{
public:
    using Distance = WeightType;

    WeightedTableDistOracle(const WeightedCRFGraph<NodeID, EdgeID, WeightType>& graph)
    {
        n = graph.num_vertices();
        const Distance inf = std::numeric_limits<Distance>::max();
        m_distances.assign(n * n, inf);

        std::vector<Distance> dist(n, inf);
        using QueueEntry = std::pair<Distance, NodeID>;
        std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> pq;

        for (NodeID source{0}; source < n; ++source)
        {
            std::fill(dist.begin(), dist.end(), inf);
            while (!pq.empty())
            {
                pq.pop();
            }

            dist[source] = static_cast<Distance>(0);
            pq.emplace(static_cast<Distance>(0), source);

            while (!pq.empty())
            {
                auto [du, u] = pq.top();
                pq.pop();
                if (du != dist[u])
                {
                    continue;
                }

                for (EdgeID e = graph.graph.vertices[u]; e < graph.graph.vertices[u + 1]; ++e)
                {
                    const NodeID v = graph.graph.edges[e];
                    const Distance w = static_cast<Distance>(graph.weights[e]);
                    if (w < static_cast<Distance>(0))
                    {
                        throw std::runtime_error("WeightedTableDistOracle requires non-negative edge weights");
                    }

                    const Distance candidate = du + w;
                    if (candidate < dist[v])
                    {
                        dist[v] = candidate;
                        pq.emplace(candidate, v);
                    }
                }
            }

            for (NodeID target{0}; target < n; ++target)
            {
                m_distances[source * n + target] = dist[target];
            }
        }
    }

    Distance get_distance(NodeID u, NodeID v) const
    {
        return m_distances[u * n + v];
    }

    void set_distance(NodeID u, NodeID v, Distance d)
    {
        m_distances[u * n + v] = d;
    }

private:
    std::vector<Distance> m_distances;
    NodeID n{0};
};
