#pragma once
#include <vector>
#include <concepts>
#include <queue>
#include "graph.hpp"

struct Noop
{
    template <typename... Args>
    void operator()(Args &&...) const noexcept {}
};

struct AlwaysFalse
{
    template <typename... Args>
    bool operator()(Args &&...) const noexcept
    {
        return false;
    }
};

template <
    typename NodeID,
    typename EdgeID,
    typename onVisited = Noop,
    typename onNotVisited,
    typename Skip = AlwaysFalse>
void bfs_single_threaded(
    const std::vector<EdgeID> &vertices,
    const std::vector<NodeID> &edges,
    const NodeID root,
    const onVisited &on_visited = {},
    const onNotVisited &on_not_visited = {},
    Skip skip = {})
    requires std::integral<NodeID> &&
             std::integral<EdgeID> &&
             std::is_invocable_v<onVisited, NodeID, NodeID> &&
             std::is_invocable_v<onNotVisited, NodeID, NodeID>
{
    std::vector<bool> visited = std::vector<bool>(vertices.size(), false);
    std::queue<NodeID> q;
    visited[root] = true;
    q.push(root);
    while (!q.empty())
    {
        NodeID current = q.front();
        q.pop();
        for (size_t i{vertices[current]}; i < vertices[current + 1]; i++)
        {
            NodeID neighbor = edges[i];
            if (skip(current, neighbor))
                continue;
            if (!visited[neighbor])
            {
                visited[neighbor] = true;
                on_not_visited(current, neighbor);
                q.push(neighbor);
            }
            else
            {
                on_visited(current, neighbor);
            }
        }
    }
}

// Generic BFS function that takes two function pointers as arguments
// to handle visited and not visited nodes.
template <typename node_T>
void bfs_parallel(
    const std::vector<size_t> &Vertices,
    const std::vector<node_T> &Edges,
    const node_T root,
    void (*const on_visited)(node_T, node_T),
    void (*const on_not_visited)(node_T, node_T))
{
    std::vector<node_T> current_frontier{Vertices.size(), -1}; // current frontier
    std::vector<node_T> next_frontier{Vertices.size(), -1};    // next frontier
    int size_current_frontier = 0;
    int size_next_frontier = 0;
    current_frontier[size_current_frontier++] = root;
    auto d = std::vector<int>(Vertices.size() - 1, -1);

#pragma omp parallel firstprivate(current_frontier, next_frontier)
    {
        std::vector<node_T> local_next_frontier{Vertices.size(), -1};
        int local_count = 0;
        while (size_current_frontier)
        {
#pragma omp for nowait
            for (int i = 0; i < size_current_frontier; i++)
            {
                int u = current_frontier[i];
                for (int j{Vertices[u]}; j < Vertices[u + 1]; j++)
                {
                    auto v = Edges[j];
                    if (d[v] < 0)
                    { // not visited
                        local_next_frontier[local_count++] = v;
                        on_not_visited(u, v);
                    }
                    else // visited
                    {
                        on_visited(u, v);
                    }
                }
            }
        }

        // calculate the offset: from what point on in the global frontier the local frontier should be copied to
        int offset;
#pragma omp critical
        {
            offset = size_next_frontier;
            size_next_frontier += local_count;
        }

        // move the local frontier into the global frontier
        for (int i{0}; i < local_count; i++)
        {
            next_frontier[i + offset] = local_next_frontier[i];
        }
        // Swap next and current frontier (there is a difference between single and master)
#pragma omp barrier
#pragma omp single
        {
            std::swap(next_frontier, current_frontier);
            size_current_frontier = size_next_frontier;
            size_next_frontier = 0;
        }
    }
}

template <typename node_T, typename edge_T>
std::vector<bool> bfs_blacklist(
    const std::vector<edge_T> &vertices,
    const std::vector<node_T> &edges,
    node_T start_node,
    std::initializer_list<std::pair<node_T, node_T>> blacklist_edges)
{
    std::vector<bool> visited{vertices.size(), false};
    std::queue<node_T> q;
    visited[start_node] = true;
    q.push(start_node);
    while (!q.empty())
    {
        node_T current = q.front();
        q.pop();
        for (size_t i{vertices[current]}; i < vertices[current + 1]; i++)
        {
            node_T neighbor = edges[i];
            if (!visited[neighbor] &&
                (std::find(blacklist_edges.begin(), blacklist_edges.end(), std::make_pair(current, neighbor)) == blacklist_edges.end()) &&
                (std::find(blacklist_edges.begin(), blacklist_edges.end(), std::make_pair(neighbor, current)) == blacklist_edges.end()))
            {
                visited[neighbor] = true;
                q.push(neighbor);
            }
        }
    }
    return visited;
}