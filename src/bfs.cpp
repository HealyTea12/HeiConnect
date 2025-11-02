#include <vector>
#include <queue>
#include "graph.hpp"

template <typename node_T>
void bfs_single_threaded(
    const std::vector<size_t> &vertices,
    const std::vector<node_T> &edges,
    const node_T root,
    void (*const on_visited)(node_T, node_T),
    void (*const on_not_visited)(node_T, node_T))
{
    std::queue<node_T> queue = {root};
    while (!queue.empty())
    {
        std::vector<bool> visited{vertices.size(), false};
        node_T current = queue.front();
        queue.pop();
        for (size_t i{vertices[current]}; i < vertices[current + 1]; ++i)
        {
            node_T neighbor = edges[i];
            if (!visited[neighbor])
            { // not visited
                on_not_visited(current, neighbor);
                queue.push(neighbor);
            }
            else
            { // visited
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

#pragma omp parallel firstprivate(current_frontier, next_frontier)
    {
        std::vector<node_T> local_next_frontier{Vertices.size(), -1};
        int local_count = 0;
        while (size_current_frontier)
        {
#pragma omp for nowait
            for (int i{0}; i < size_current_frontier; i++)
            {
                int u = current_frontier[i];
                for (int j{Vertices[u]}; j < Vertices[u + 1]; j++)
                {
                    v = Edges[j];
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
std::unordered_set<node_T> bfs_ignore_edges(
    const std::vector<edge_T> &Vertices,
    const std::vector<node_T> &Edges,
    const node_T root,
    const std::pair<node_T, node_T> &ignore_edges...)
{
    std::unordered_set<node_T> visited;
    bfs_single_threaded<node_T>(
        Vertices,
        Edges,
        root,
        [&](node_T u, node_T v)
        {
            // on visited
        },
        [&](node_T u, node_T v)
        {
            // on not visited
            if (!((u == ignore_edges.first && v == ignore_edges.second) ||
                  (u == ignore_edges.second && v == ignore_edges.first)))
            {
                visited.insert(v);
            }
        });
    return visited;
}