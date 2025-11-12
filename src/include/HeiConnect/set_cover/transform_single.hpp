#pragma once
#include <concepts>
#include <vector>
#include <queue>
#include <unordered_set>
#include <cassert>

#include "min_cut/simple_mincut.hpp"
#include "set_cover/set_cover.hpp"
#include "bfs.hpp"

// node_T/edge_T is a generic type that indexes nodes/edges e.g unsigned int
template <typename node_T>
struct Edge
{
    node_T u;
    node_T v;
    double weight;
};

// template <typename G>
// concept graph = requires(G &graph) {
//     graph.num_nodes()->std::integral_type
// };

struct pair_hash
{
    inline std::size_t operator()(const std::pair<int, int> &v) const
    {
        return v.first * 31 + v.second;
    }
};

// returns the edge idx between u and v.
// Assumes that there is an edge between u and v.
// If there is no such edge, assert false.
// Could consider adding an optional return type instead.
// But, in this context, if the edge doesn't exist, then there is a bug in the code.
template <typename NodeID, typename EdgeID>
static inline EdgeID get_edge_index(
    const std::vector<EdgeID> &vertices,
    const std::vector<NodeID> &edges,
    NodeID u,
    NodeID v) noexcept
{
    for (auto i{vertices[u]}; i < vertices[u + 1]; i++)
    {
        if (edges[i] == v)
        {
            return i;
        }
    }
    assert(false);
}

template <typename node_T, typename edge_T>
    requires std::integral<node_T> && std::integral<edge_T>
SetCover construct_set_cover(
    const std::vector<edge_T> &vertices,
    const std::vector<node_T> &edges,
    const std::vector<double> &weights,
    const std::vector<size_t> &link_vertices,
    const std::vector<size_t> &link_edges,
    const std::vector<double> &link_weights)
{
    // auto cycle_edges = std::unordered_set<std::pair<node_T, node_T>>{0, pair_hash};
    // for checking during the second BFS
    auto cycle_edges = std::vector<bool>(edges.size(), false);
    std::vector<std::pair<node_T, node_T>> cycle_edge_vec{};
    std::vector<int> distances = std::vector<int>(vertices.size(), 0);
    std::vector<bool> visited_nodes = std::vector<bool>(vertices.size(), false);
    std::vector<node_T> parent = std::vector<node_T>(vertices.size(), 0);
    // BFS to find cycles and rooted tree
    std::queue<std::pair<node_T, int>> q;
    node_T root = static_cast<node_T>(0);
    parent[root] = root;
    q.push({root, 0});
    visited_nodes[root] = true;
    while (!q.empty())
    {
        auto curr = q.front();
        node_T current_node = curr.first;
        auto depth = curr.second;
        q.pop();
        for (size_t i{vertices[current_node]}; i < vertices[current_node + 1]; i++)
        {
            node_T neighbor = edges[i];
            if (!visited_nodes[neighbor])
            {
                visited_nodes[neighbor] = true;
                parent[neighbor] = current_node;
                distances[neighbor] = depth + 1;
                q.push({neighbor, depth + 1});
            }
            else
            {
                if (distances[neighbor] < depth)
                    continue; // skip back edges to ancestors
                cycle_edge_vec.emplace_back(current_node, neighbor);
            }
        }
    }
    // for each cycle edge, reconstruct the cycle
    // could rethink data struct
    auto cycles = std::vector<std::vector<std::pair<node_T, node_T>>>(cycle_edge_vec.size());
    assert(cycle_edge_vec.size() == cycles.size());
    for (size_t i{}; i < cycle_edge_vec.size(); i++)
    {
        auto &cycle = cycles[i];
        cycle.emplace_back(cycle_edge_vec[i]);
        node_T u = cycle_edge_vec[i].first;
        node_T v = cycle_edge_vec[i].second;
        cycle_edges[get_edge_index(vertices, edges, u, v)] = true;
        cycle_edges[get_edge_index(vertices, edges, v, u)] = true;
        if (distances[u] < distances[v])
        {
            cycle_edges[get_edge_index(vertices, edges, v, parent[v])] = true;
            cycle_edges[get_edge_index(vertices, edges, parent[v], v)] = true;
            v = parent[v];
        }
        while (u != v)
        {
            cycle.emplace_back(u, parent[u]);
            cycle.emplace_back(v, parent[v]);
            // could store the indices, so you don't have to look them up
            cycle_edges[get_edge_index(vertices, edges, u, parent[u])] = true;
            cycle_edges[get_edge_index(vertices, edges, parent[u], u)] = true;
            cycle_edges[get_edge_index(vertices, edges, v, parent[v])] = true;
            cycle_edges[get_edge_index(vertices, edges, parent[v], v)] = true;
            u = parent[u];
            v = parent[v];
        }
    }

    // iterate over all edges to find min cuts
    auto min_cut = global_mincut_simple({{vertices, edges},
                                         weights}); // should substitute with cactus min cut
    std::vector<bool> min_cuts{};
    size_t n_min_cuts{};
    size_t n_vertices = vertices.size() - 1;
    for (node_T u{}; u < vertices.size() - 1; u++)
    {
        for (edge_T i{vertices[u]}; i < vertices[u + 1]; i++)
        {
            node_T v = edges[i];
            if (u < v && !cycle_edges[i])
            {
                if (weights[i] == min_cut)
                {
                    min_cuts.resize((n_min_cuts + 1) * n_vertices, false);
                    min_cuts[n_min_cuts * n_vertices + u] = true;
                    bfs_single_threaded(
                        vertices,
                        edges,
                        u,
                        [](node_T from, node_T to) noexcept {},
                        [&min_cuts, n_min_cuts, n_vertices](node_T from, node_T to) noexcept
                        {
                            min_cuts[n_min_cuts * n_vertices + to] = true;
                        },
                        [u, v](node_T from, node_T to) noexcept
                        {
                            if (from == u && to == v)
                                return true;
                            return false;
                        });
                    n_min_cuts++;
                }
            }
        }
    }
    for (const auto &cycle : cycles)
    {
        for (size_t i{}; i < cycle.size() - 1; i++)
        {
            for (size_t j = i + 1; j < cycle.size(); j++)
            {
                auto edge1 = cycle[i];
                auto edge2 = cycle[j];
                min_cuts.resize((n_min_cuts + 1) * vertices.size(), false);
                min_cuts[n_min_cuts * n_vertices + edge1.first] = true;
                bfs_single_threaded(
                    vertices,
                    edges,
                    edge1.first,
                    [](node_T from, node_T to) noexcept {},
                    [&min_cuts, n_min_cuts, n_vertices](node_T from, node_T to) noexcept
                    {
                        min_cuts[n_min_cuts * n_vertices + to] = true;
                    },
                    [&edge1, &edge2](node_T from, node_T to) noexcept
                    {
                        if ((from == edge1.first && to == edge1.second) ||
                            (from == edge2.first && to == edge2.second) ||
                            (from == edge2.second && to == edge2.first))
                            return true;
                        return false;
                    });
            }
        }
    }
    // for each link, determine which cuts it crosses
    std::vector<size_t> a = std::vector<size_t>(link_edges.size() + 1, static_cast<size_t>(0));
    std::vector<size_t> b{};
    for (size_t i{}; i < link_vertices.size(); i++)
    {
        for (size_t j{link_vertices[i]}; j < link_vertices[i + 1]; j++)
        {
            auto link = Edge{i, link_edges[j], link_weights[j]};
            for (size_t k{}; k < n_min_cuts; k++)
            {
                if (min_cuts[k * n_vertices + link.u] != min_cuts[k * n_vertices + link.v])
                {
                    b.push_back(k);
                }
            }
            a[j + 1] = b.size();
        }
    }
    return SetCover{a, b, link_weights};
}