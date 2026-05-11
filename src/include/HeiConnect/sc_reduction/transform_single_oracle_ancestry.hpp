#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <limits>
#include <stack>
#include <utility>
#include <vector>

#include "HeiConnect/sc_reduction/cactus_min_cuts.hpp"
#include "HeiConnect/sc_reduction/transform_single_partition_matrix_utils.hpp"

using ull = unsigned long long;
using timer_type = uint32_t;

namespace HeiConnect_details
{
    template <class node_T>
        requires std::unsigned_integral<node_T>
    struct TinToutResult
    {
        /// @brief when explored in DFS
        std::vector<timer_type> tin;
        /// @brief when finished exploring in DFS (when it pops)
        std::vector<timer_type> tout;
        // currently from last explored node to first
        std::vector<std::vector<node_T>> cycles;
        /// @brief =1 if edge is a cycle edge, 0 otherwise
        std::vector<char> is_cycle_edge;
        std::vector<node_T> parent;
    };

    template <typename node_T, typename edge_T>
        requires std::integral<node_T> && std::integral<edge_T>
    TinToutResult<node_T>
    dfs_tin_tout_cycles(
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges)
    {
        // maybe having a vector of structs is faster
        auto tin = std::vector<timer_type>(vertices.size() - 1, -1);
        auto tout = std::vector<timer_type>(vertices.size() - 1, -1);
        auto cycle_edges = std::vector<std::pair<node_T, node_T>>{};
        auto parent = std::vector<node_T>(vertices.size() - 1, -1);
        auto is_cycle_edge = std::vector<char>(edges.size(), 0);
        // 0 = unseen, 1 = in current DFS stack, 2 = fully processed
        auto state = std::vector<char>(vertices.size() - 1, 0);
        auto next_edge = std::vector<edge_T>(vertices.size() - 1, 0);
        auto stack = std::stack<node_T>{};
        node_T root = static_cast<node_T>(0); // can be any node
        timer_type timer{0};
        stack.push(root);
        parent[root] = root;
        while (!stack.empty())
        {
            node_T current_node = stack.top();

            if (state[current_node] == 0)
            {
                tin[current_node] = timer++;
                state[current_node] = 1;
                next_edge[current_node] = vertices[current_node];
            }

            bool advanced = false;
            for (edge_T &i = next_edge[current_node]; i < vertices[current_node + 1]; ++i)
            {
                node_T neighbor = edges[i];
                if (neighbor == parent[current_node])
                {
                    continue;
                }

                if (state[neighbor] == 0)
                {
                    parent[neighbor] = current_node;
                    stack.push(neighbor);
                    ++i;
                    advanced = true;
                    break;
                }

                // In undirected DFS, only back-edges to active ancestors represent cycles.
                if (state[neighbor] == 1 && tin[neighbor] < tin[current_node])
                {
                    cycle_edges.emplace_back(current_node, neighbor);
                }
            }

            if (!advanced)
            {
                tout[current_node] = timer++;
                state[current_node] = 2;
                stack.pop();
            }
        }

        std::vector<std::vector<node_T>> cycles{};
        for (const auto &edge : cycle_edges)
        {
            std::vector<node_T> cycle{};
            node_T u = edge.first;
            node_T v = edge.second;
            // mark cycle edges
            // TODO: optimize by storing indices during DFS
            is_cycle_edge[get_edge_index(vertices, edges, u, v)] = 1;
            is_cycle_edge[get_edge_index(vertices, edges, v, u)] = 1;
            while (u != v)
            {
                // mark cycle edges
                // TODO: optimize by storing indices during DFS
                is_cycle_edge[get_edge_index(vertices, edges, u, parent[u])] = 1;
                is_cycle_edge[get_edge_index(vertices, edges, parent[u], u)] = 1;
                cycle.emplace_back(u);
                u = parent[u];
            }
            cycle.emplace_back(v);
            cycles.emplace_back(cycle);
        }

        return {tin, tout, cycles, is_cycle_edge, parent};
    }

    template <typename node_T, typename edge_T, typename weight_T>
        requires std::integral<node_T> && std::integral<edge_T>
    weight_T calculate_cactus_min_cut(
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges,
        const std::vector<weight_T> &weights,
        const std::vector<std::vector<node_T>> &cycles,
        const std::vector<char> &is_cycle_edge)
    {
        weight_T min_cut = std::numeric_limits<weight_T>::max();
        // iterate over all edges to find min cuts
        for (node_T u{}; u < vertices.size() - 1; u++)
        {
            for (edge_T e{vertices[u]}; e < vertices[u + 1]; e++)
            {
                if (!is_cycle_edge[e] && weights[e] < min_cut)
                {
                    min_cut = weights[e];
                }
            }
        }
        for (const auto &cycle : cycles)
        {
            for (auto i = 0; i < cycle.size() - 1; i++)
            {
                node_T edge1_u = cycle[i];
                node_T edge1_v = cycle[(i + 1) % cycle.size()];
                auto edge1_idx = get_edge_index(vertices, edges, edge1_u, edge1_v);
                for (auto j = i + 1; j < cycle.size(); j++)
                {
                    node_T edge2_u = cycle[j];
                    node_T edge2_v = cycle[(j + 1) % cycle.size()];
                    auto edge2_idx = get_edge_index(vertices, edges, edge2_u, edge2_v);
                    weight_T cycle_cut = weights[edge1_idx] + weights[edge2_idx];
                    if (cycle_cut < min_cut)
                    {
                        min_cut = cycle_cut;
                    }
                }
            }
        }
        return min_cut;
    }

    template <typename node_T, typename edge_T, typename weight_T>
        requires std::integral<node_T> && std::integral<edge_T>
    std::vector<node_T> calculate_tree_cactus_min_cuts(
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges,
        const std::vector<weight_T> &weights,
        const std::vector<char> &is_cycle_edge,
        weight_T min_cut,
        std::vector<node_T> &parent,
        node_T root)
    {
        std::vector<node_T> tree_mcs{};
        for (node_T u = 0; u < vertices.size() - 1; u++)
        {
            for (edge_T e{vertices[u]}; e < vertices[u + 1]; e++)
            {
                node_T v = edges[e];
                if (!is_cycle_edge[e] && weights[e] == min_cut)
                {
                    if (u == parent[v])
                    {
                        tree_mcs.emplace_back(v);
                    }
                    else
                    {
                        tree_mcs.emplace_back(u);
                    }
                }
            }
        }
        return tree_mcs;
    }

    template <typename node_T, typename edge_T, typename weight_T>
        requires std::integral<node_T> && std::integral<edge_T>
    CactusMinCuts<node_T> calculate_all_cactus_min_cuts(
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges,
        const std::vector<weight_T> &weights,
        const std::vector<std::vector<node_T>> &cycles,
        const std::vector<char> &is_cycle_edge,
        weight_T min_cut,
        std::vector<node_T> &parent,
        node_T root)

    {
        // TODO: might be adding duplicate cuts
        // don't forget to squeeze
        std::vector<node_T> tree_mcs = calculate_tree_cactus_min_cuts(
            vertices, edges, weights, is_cycle_edge, min_cut, parent, root);

        // There is an important exception case to take care of here
        // We assume that the cycles are stored in traversal order
        // Each cycle mincut is the child of the first and second edges
        std::vector<std::pair<node_T, node_T>> cycle_mcs{};
        // TODO: THIS TEMPORARILY REVERSES THE CYCLE ORDER FOR DEBUGGING
        auto reversed_cycles = cycles;
        for (auto &cycle : reversed_cycles)
        {
            std::reverse(cycle.begin(), cycle.end());
        }
        for (const auto &cycle : reversed_cycles)
        {
            assert(cycle.size() >= 3 && "Cycles must have at least 3 nodes; bug in cycle detection");
            for (size_t i = 0; i < cycle.size() - 2; i++)
            {
                node_T edge1_u = cycle[i];
                node_T edge1_v = cycle[i + 1];
                auto edge1_idx = get_edge_index(vertices, edges, edge1_u, edge1_v);
                for (size_t j = i + 1; j < cycle.size() - 1; j++)
                {
                    node_T edge2_u = cycle[j];
                    node_T edge2_v = cycle[j + 1];
                    auto edge2_idx = get_edge_index(vertices, edges, edge2_u, edge2_v);
                    weight_T cycle_cut = weights[edge1_idx] + weights[edge2_idx];
                    if (cycle_cut == min_cut)
                    {
                        cycle_mcs.emplace_back(edge1_v, edge2_v);
                    }
                }
            }
            // the last edge of the cycle is an exception as it is not part of the tree
            node_T f_u{cycle[cycle.size() - 1]};
            node_T f_v{cycle[0]};
            node_T f_e = get_edge_index(vertices, edges, f_u, f_v);
            weight_T f_w = weights[f_e];
            for (size_t i{0}; i < cycle.size() - 1; i++)
            {
                node_T u = cycle[i];
                node_T v = cycle[i + 1];
                edge_T e = get_edge_index(vertices, edges, u, v);
                weight_T w = weights[e];
                if (w + f_w == min_cut)
                {
                    tree_mcs.emplace_back(v);
                }
            }
        }
        return CactusMinCuts{tree_mcs, cycle_mcs};
    }

    // is u ancestor of v
    template <typename node_T>
        requires std::integral<node_T>
    bool isAncestor(
        const std::vector<timer_type> &tin,
        const std::vector<timer_type> &tout,
        node_T u,
        node_T v)
    {
        return tin[u] <= tin[v] && tout[u] >= tout[v];
    }

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
