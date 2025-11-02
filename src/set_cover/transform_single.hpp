#pragma once
#include <concepts>
#include <vector>
#include <queue>
#include <unordered_set>
#include "set_cover/set_cover.cpp"

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
    std::unordered_set<std::pair<node_T, node_T>> cycle_edges;
    std::vector<std::pair<node_T, node_T>> cycle_edge_vec;
    std::vector<int> distances{vertices.size(), 0};
    std::vector<bool> visited_nodes{vertices.size(), false};
    std::vector<node_T> parent{vertices.size(), 0};
    std::queue<std::pair<node_T, int>> q;
    node_T root = static_cast<node_T>(0);
    parent[root] = root;
    q.push({root, 0});
    while (!q.empty())
    {
        auto [current_node, depth] = q.front();
        q.pop();
        for (size_t i{vertices[current_node]}; i < vertices[current_node + 1]; i++)
        {
            node_T neighbor = edges[i];
            if (!visited_nodes[neighbor])
            {
                visited_nodes[neighbor] = true;
                parent[neighbor] = current_node;
                q.push({neighbor, depth + 1});
            }
            else
            {
                cycle_edge_vec.emplace_back({current_node, neighbor});
            }
        }
    }
    std::vector<std::vector<std::pair<node_T, node_T>>> cycles{cycle_edge_vec.size()};
    for (size_t i{}; i < cycle_edge_vec.size(); i++)
    {
        auto cycle = cycles[i];
        cycle.emplace_back(cycle_edge_vec[i]);
        cycle_edges.insert(cycle_edge_vec[i]);
        node_T u = cycle_edge_vec[i].first;
        node_T v = cycle_edge_vec[i].second;
        if (distances[u] < distances[v])
        {
            cycle.emplace_back(v, parent[v]);
            v = parent[v];
        }
        while (u != v)
        {
            cycle.emplace_back(u, parent[u]);
            cycle.emplace_back(v, parent[v]);
            cycle_edges.insert({u, parent[u]});
            cycle_edges.insert({v, parent[v]});
            u = parent[u];
            v = parent[v];
        }
    }

    // iterate over all edges to find min cuts
    std::vector<std::vector<bool>> min_cuts;
    for (node_T u{}; u < vertices.size(); u++)
    {
        for (auto i{vertices[u]}; i < vertices[u + 1]; i++)
        {
            node_T v = edges[i];
            if (cycle_edges.find({u, v}) == cycle_edges.end() &&
                cycle_edges.find({v, u}) == cycle_edges.end())
            {
                if (weights[i] == min_cut)
                {
                    min_cuts.push_back(bfs_blacklist(
                        vertices,
                        edges,
                        u,
                        {u, v}))
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
                min_cuts.push_back(bfs_blacklist(
                    vertices,
                    edges,
                    edge1.first,
                    {edge1, edge2}));
            }
        }
    }
    // for each link, determine which cuts it crosses
    std::vector<int> a{links.size(), 0};
    std::vector<int> b{};
    for (size_t i{}; i < links.size(); i++)
    {
        auto link = links[i];
        for (size_t j{}; j < min_cuts.size(); j++)
        {
            auto cut = min_cuts[j];
            if (cut[link.u] != cut[link.v])
            {
                b.push_back(j);
            }
        }
        a[i] = b.size();
    }
    return SetCover{a, b, link_weights};
}