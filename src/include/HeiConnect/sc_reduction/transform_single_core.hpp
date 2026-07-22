#pragma once
#include <concepts>
#include <vector>
#include <queue>
#include <tuple>
#include <stack>
#include <unordered_set>
#include <cassert>
#include <limits>
#include <omp.h>
#include <immintrin.h>
#include <type_traits>

#include "HeiConnect/min_cut/simple_mincut.hpp"
#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/bfs.hpp"
#include "HeiConnect/sc_reduction/cactus_min_cuts.hpp"

using ull = unsigned long long;

namespace HeiConnect_details
{
    template <typename node_T, typename edge_T>
        requires std::unsigned_integral<node_T> && std::unsigned_integral<edge_T>
    struct Info
    {
        node_T parent;
        edge_T parent_edge_index; // from parent to v
        int distance;
    };

    template <typename node_T, typename edge_T>
        requires std::integral<node_T> && std::integral<edge_T>
    static std::tuple<
        std::vector<std::pair<node_T, node_T>>,
        std::vector<Info<node_T, edge_T>>>
    find_cycles_and_root_tree(
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges,
        const std::vector<double> &weights,
        const std::vector<size_t> &link_vertices,
        const std::vector<size_t> &link_edges,
        const std::vector<double> &link_weights)
    {
        std::vector<Info<node_T, edge_T>> info = std::vector<Info<node_T, edge_T>>(vertices.size() - 1);
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
                    if (info[neighbor].distance == depth && neighbor < current_node)
                        continue; // count same-depth undirected edges only once
                    cycle_edge_vec.emplace_back(current_node, neighbor);
                }
            }
        }
        return std::tuple<
            std::vector<std::pair<node_T, node_T>>,
            std::vector<Info<node_T, edge_T>>>{cycle_edge_vec, info};
    }
};
