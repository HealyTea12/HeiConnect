#pragma once
#include <concepts>
#include <vector>
#include <queue>
#include <unordered_set>
#include <cassert>
#include <omp.h>

#include "HeiConnect/min_cut/simple_mincut.hpp"
#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/bfs.hpp"

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
    return std::numeric_limits<EdgeID>::max(); // max index is reserved for invalid
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

    double start = omp_get_wtime();
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
    double end = omp_get_wtime();
    std::cout << "BFS to find cycles took " << (end - start) << " seconds." << std::endl;
    // for each cycle edge, reconstruct the cycle
    // could rethink data struct
    start = omp_get_wtime();
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
    end = omp_get_wtime();
    std::cout << "Cycle reconstruction took " << (end - start) << " seconds." << std::endl;
    // iterate over all edges to find min cuts
    // auto min_cut = global_mincut_simple({{vertices, edges},
    //                                     weights}); // should substitute with cactus min cut
    // this part is awful
    start = omp_get_wtime();
    auto min_cut = std::numeric_limits<double>::max();
    for (node_T u{}; u < vertices.size() - 1; u++)
    {
        for (edge_T e{vertices[u]}; e < vertices[u + 1]; e++)
        {
            if (!cycle_edges[e] && weights[e] < min_cut)
            {
                min_cut = weights[e];
            }
        }
    }
    for (const auto &cycle : cycles)
    {
        for (auto i = 0; i < cycle.size() - 1; i++)
        {
            auto edge1 = cycle[i];
            auto edge1_idx = get_edge_index(vertices, edges, edge1.first, edge1.second);
            for (auto j = i + 1; j < cycle.size(); j++)
            {
                auto edge2 = cycle[j];
                auto edge_idx = get_edge_index(vertices, edges, edge2.first, edge2.second);
                double cycle_cut = weights[edge1_idx] + weights[edge_idx];
                if (cycle_cut < min_cut)
                {
                    min_cut = cycle_cut;
                }
            }
        }
    }
    end = omp_get_wtime();
    std::cout << "Cactus min cut computation took " << (end - start) << " seconds." << std::endl;

    // find and partition min cuts
    start = omp_get_wtime();
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
    end = omp_get_wtime();
    std::cout << "Min cut tree edge partitioning took " << (end - start) << " seconds." << std::endl;
    start = omp_get_wtime();
    for (const auto &cycle : cycles)
    {
        for (size_t i{}; i < cycle.size() - 1; i++)
        {
            for (size_t j = i + 1; j < cycle.size(); j++)
            {
                auto edge1 = cycle[i];
                auto edge2 = cycle[j];
                auto edge1_idx = get_edge_index(vertices, edges, edge1.first, edge1.second);
                auto edge2_idx = get_edge_index(vertices, edges, edge2.first, edge2.second);
                double cycle_cut = weights[edge1_idx] + weights[edge2_idx];
                if (cycle_cut == min_cut)
                {
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
                    n_min_cuts++;
                }
            }
        }
    }
    end = omp_get_wtime();
    std::cout << "min cut cycle edge partitioning took " << (end - start) << " seconds." << std::endl;
    // for each link, determine which cuts it crosses
    start = omp_get_wtime();
    std::vector<size_t> a = std::vector<size_t>(link_edges.size() + 1, static_cast<size_t>(0));
    std::vector<size_t> b{};
    b.reserve(link_edges.size() * n_min_cuts);
    for (size_t u{}; u < link_vertices.size() - 1; u++)
    {
        for (size_t e{link_vertices[u]}; e < link_vertices[u + 1]; e++)
        {
            // auto link = Edge{u, link_edges[e], link_weights[e]};
            node_T v = link_edges[e];
            for (size_t k{}; k < n_min_cuts; k++)
            {
                if (min_cuts[k * n_vertices + u] != min_cuts[k * n_vertices + v])
                {
                    b.emplace_back(k);
                }
            }
        }
        a[u + 1] = b.size();
    }
    end = omp_get_wtime();
    std::cout << "Link cut determination/setcover construction took " << (end - start) << " seconds." << std::endl;
    // this is probably copying, we should change to move
    return SetCover{a, b, link_weights};
}

template <typename node_T, typename edge_T>
    requires std::integral<node_T> && std::integral<edge_T>
static auto find_cycles_and_root_tree(
    const std::vector<edge_T> &vertices,
    const std::vector<node_T> &edges,
    const std::vector<double> &weights,
    const std::vector<size_t> &link_vertices,
    const std::vector<size_t> &link_edges,
    const std::vector<double> &link_weights)
{
    struct Info
    {
        node_T parent;
        edge_T parent_edge_index; // from parent to v
        int distance;
    };
    std::vector<Info> info = std::vector<Info>(vertices.size() - 1);
    std::vector<bool> visited_nodes = std::vector<bool>(vertices.size() - 1, false);
    // store parent as (parent_node, edge_index) edge_index is for the edge going from v to parent(v)
    // otherwise, we would have to do a linear search over the neighbours to find the edge index
    // when reconstructing the cycles
    std::queue<std::pair<node_T, int>> q;
    node_T root = static_cast<node_T>(0);
    info[root] = {root, std::numeric_limits<edge_T>::max(), 0};
    visited_nodes[root] = true;
    // parent, child
    std::vector<std::pair<node_T, node_T>> cycle_edge_vec{};
    cycle_edge_vec.reserve(edges.size());
    q.push({root, 0});
    while (!q.empty())
    {
        auto curr = q.front();
        node_T current_node = curr.first;
        auto depth = curr.second;
        q.pop();
        for (edge_T i{vertices[current_node]}; i < vertices[current_node + 1]; i++)
        {
            node_T neighbor = edges[i];
            if (!visited_nodes[neighbor])
            {
                visited_nodes[neighbor] = true;
                info[neighbor] = {current_node, i, depth + 1};
                q.push({neighbor, depth + 1});
            }
            else
            {
                if (info[neighbor].distance < depth)
                    continue; // skip back edges to ancestors
                cycle_edge_vec.emplace_back(current_node, neighbor);
            }
        }
    }
    return std::tuple<
        std::vector<std::pair<node_T, node_T>>,
        std::vector<Info>>{cycle_edge_vec, info};
}

// template <typename node_T, typename edge_T>
//     requires std::integral<node_T> && std::integral<edge_T>
// SetCover construct_set_cover1(
//     const std::vector<edge_T> &vertices,
//     const std::vector<node_T> &edges,
//     const std::vector<double> &weights,
//     const std::vector<size_t> &link_vertices,
//     const std::vector<size_t> &link_edges,
//     const std::vector<double> &link_weights)
// {
//     // auto cycle_edges = std::unordered_set<std::pair<node_T, node_T>>{0, pair_hash};
//     // for checking during the second BFS
//     auto cycle_edges = std::vector<bool>(edges.size(), false);
//     // BFS to find cycles and rooted tree
//     std::vector<std::pair<node_T, node_T>> cycle_edge_vec{};
//     std::vector<find_cycles_and_root_tree::Info> info{};
//     auto [cycle_edge_vec, info] = find_cycles_and_root_tree(vertices, edges, weights, link_vertices, link_edges, link_weights);
//     // for each cycle edge, reconstruct the cycle
//     // could rethink data struct
//     auto cycles = std::vector<std::vector<std::pair<node_T, node_T>>>(cycle_edge_vec.size());
//     for (size_t i{}; i < cycle_edge_vec.size(); i++)
//     {
//         auto &cycle = cycles[i];
//         cycle.emplace_back(cycle_edge_vec[i]);
//         node_T u = cycle_edge_vec[i].first;
//         node_T v = cycle_edge_vec[i].second;
//         // mark only one direction, as we only traverse edges in this direction during min cut
//         cycle_edges[info[v].parent_edge_index] = true;
//         if (info[u].distance < info[v].distance)
//         {
//             cycle_edges[info[v].parent_edge_index] = true;
//             v = info[v].parent_node;
//         }
//         while (u != v)
//         {
//             cycle.emplace_back(info[u].parent_node, u);
//             cycle.emplace_back(info[v].parent_node, v);
//             cycle_edges[info[u].parent_edge_index] = true;
//             cycle_edges[info[v].parent_edge_index] = true;
//             u = info[u].parent_node;
//             v = info[v].parent_node;
//         }
//     }
//
//     // iterate over all edges to find min cuts
//     // auto min_cut = global_mincut_simple({{vertices, edges},
//     //                                     weights}); // should substitute with cactus min cut
//     // this part is awful
//     auto min_cut = std::numeric_limits<double>::max();
//     bfs_single_threaded(
//         vertices,
//         edges,
//         static_cast<node_T>(0),
//         [&cycle_edges, &info, &weights, &min_cut](node_T from, node_T to) noexcept
//         {
//             if (!cycle_edges[info[to].parent_edge_index] && weights[info[to].parent_edge_index] < min_cut)
//             {
//                 min_cut = weights[info[to].parent_edge_index];
//             }
//         },
//         [&cycle_edges, &info, &weights, &min_cut](node_T from, node_T to) noexcept
//         {
//             if (!cycle_edges[info[to].parent_edge_index] && weights[info[to].parent_edge_index] < min_cut)
//             {
//                 min_cut = weights[info[to].parent_edge_index];
//             }
//         },
//         [](node_T from, node_T to) noexcept
//         {
//             return false;
//         });
//     for (const auto &cycle : cycles)
//     {
//         for (auto i = 0; i < cycle.size() - 1; i++)
//         {
//             auto edge1 = cycle[i];
//             auto edge1_idx = info[edge1.second].parent_edge_index;
//             for (auto j = i + 1; j < cycle.size(); j++)
//             {
//                 auto edge2 = cycle[j];
//                 auto edge_idx = info[edge2.second].parent_edge_index;
//                 double cycle_cut = weights[edge1_idx] + weights[edge_idx];
//                 if (cycle_cut < min_cut)
//                 {
//                     min_cut = cycle_cut;
//                 }
//             }
//         }
//     }
//
//     // wrong
//     std::vector<uint8_t> min_cuts{};
//     size_t n_min_cuts{};
//     size_t n_vertices = vertices.size() - 1;
//     min_cuts.resize(n_vertices * n_vertices * (n_vertices - 1) / (2 * 8), 0);
//     for (node_T u{}; u < vertices.size() - 1; u++)
//     {
//         for (edge_T i{vertices[u]}; i < vertices[u + 1]; i++)
//         {
//             node_T v = edges[i];
//             if (u < v && !cycle_edges[i])
//             {
//                 if (weights[i] == min_cut)
//                 {
//                     auto i = n_min_cuts * n_vertices + u;
//                     min_cuts[i / 8] &= 1 << (i % 8);
//                     bfs_single_threaded(
//                         vertices,
//                         edges,
//                         u,
//                         [](node_T from, node_T to) noexcept {},
//                         [&min_cuts, n_min_cuts, n_vertices](node_T from, node_T to) noexcept
//                         {
//                             auto i = n_min_cuts * n_vertices + to;
//                             min_cuts[i / 8] &= 1 << (i % 8);
//                         },
//                         [u, v](node_T from, node_T to) noexcept
//                         {
//                             if (from == u && to == v)
//                                 return true;
//                             return false;
//                         });
//                     n_min_cuts++;
//                 }
//             }
//         }
//     }
//     for (const auto &cycle : cycles)
//     {
//         for (size_t i{}; i < cycle.size() - 1; i++)
//         {
//             for (size_t j = i + 1; j < cycle.size(); j++)
//             {
//                 auto edge1 = cycle[i];
//                 auto edge2 = cycle[j];
//                 auto edge1_idx = get_edge_index(vertices, edges, edge1.first, edge1.second);
//                 auto edge2_idx = get_edge_index(vertices, edges, edge2.first, edge2.second);
//                 double cycle_cut = weights[edge1_idx] + weights[edge2_idx];
//                 if (cycle_cut == min_cut)
//                 {
//                     auto i = n_min_cuts * n_vertices + edge1.first;
//                     min_cuts[i / 8] &= 1 << (i % 8);
//                     bfs_single_threaded(
//                         vertices,
//                         edges,
//                         edge1.first,
//                         [](node_T from, node_T to) noexcept {},
//                         [&min_cuts, n_min_cuts, n_vertices](node_T from, node_T to) noexcept
//                         {
//                             auto i = n_min_cuts * n_vertices + to;
//                             min_cuts[i / 8] &= 1 << (i % 8);
//                         },
//                         [&edge1, &edge2](node_T from, node_T to) noexcept
//                         {
//                             if ((from == edge1.first && to == edge1.second) ||
//                                 (from == edge2.first && to == edge2.second) ||
//                                 (from == edge2.second && to == edge2.first))
//                                 return true;
//                             return false;
//                         });
//                     n_min_cuts++;
//                 }
//             }
//         }
//     }
//     // for each link, determine which cuts it crosses
//     std::vector<size_t> a = std::vector<size_t>(link_edges.size() + 1, static_cast<size_t>(0));
//     std::vector<size_t> b{};
//     b.reserve(link_edges.size() * 3);
//     for (size_t i{}; i < link_vertices.size() - 1; i++)
//     {
//         for (size_t j{link_vertices[i]}; j < link_vertices[i + 1]; j++)
//         {
//             auto link = Edge{i, link_edges[j], link_weights[j]};
//             for (size_t k{}; k < n_min_cuts; k++)
//             {
//                 auto u_idx = k * n_vertices + link.u;
//                 auto v_idx = k * n_vertices + link.v;
//                 auto cond = (min_cuts[u_idx / 8] & (1 << (u_idx % 8))) ^ (min_cuts[v_idx / 8] & (1 << (v_idx % 8)));
//                 if (cond)
//                 {
//                     b.emplace_back(k);
//                 }
//             }
//             a[j + 1] = b.size();
//         }
//     }
//     return SetCover{a, b, link_weights};
// }