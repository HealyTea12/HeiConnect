#pragma once

#include <algorithm>
#include <chrono>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include "HeiConnect/conn_aug/reducers/common.hpp"
#include "HeiConnect/data_structures/distance_oracle.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/pipeline/common.hpp"
#include "HeiConnect/tools/timer.hpp"


/* TODO:
    1. Fix tree_path so that it ignores the cycle nodes in the cactus tree of the cactus graph.
    2. When projecting in, project in will try to set_distance to a cycle node, which doesn't exist in the distance
   oracle.
   3. In projectOut, we need to make sure that the distance between nodes in the original graph skip the cycle nodes
   (maybe have the graph be weighted and give the edges connected to the cycle node cost 0?)
   4. projectIn is incorrect. We need to process paths in order of length and project to immediate subintervals.

   (low priority):
   1. distance sum might overflow. This won't be a problem for small distances, but keep in mind
 */


/*
  The goal of this reducer is to find all many to one link reductions in the connectivity augmentation problem.
  A set of links L_dom dominates another link (u,v) if L_dom covers the path from u to v in the original graph
  (cutting through any cycles) and w(L_dom) <= w(u,v).

  We assume that the input graph is already a cactus graph.

  1. Find all shortest paths for all pairs of nodes in the link graph.
  Do this with Dijkstra's or Floyd-Warshall.

  2. Propagate the shortest path distances "inwards" in the original graph.
  Consider a path u_0 -> u_1 -> ... -> u_k in the original graph.
  Then every sub path P = u_i -> ... -> u_j, i < j, will get the distance min_{x <= i, y >= j} d(u_x, u_y).
  In other words, if there exists a path that encapsulates the sub path P
  and has a shorter distance in the link graph, then we can assign that distance to P.

  To achieve this we propagate the distances from the leaves of the original graph inwards towards the root.

  3. Propagate the shortest path distances "outwards" in the original graph. (need a better name for this step)
    Consider a path P = u_0 -> u_1 -> ... -> u_k in the original graph.
    If there exists a combination of paths P_1, P_2, ..., P_m that cover P and sum(d(P_i)) <= d(P),
    then we can assign d(P) = sum(d(P_i)).

    For this we start with paths of length 1 and propagate outwards to their neighbouring paths.
*/
template<int RecordStatsLevel = 0>
class FullReducer
{
public:
    static std::string name()
    {
        return "Full Connectivity Augmentation Reducer";
    }

    FullReducer(bool project_in, bool project_out) : m_projectIn(project_in), m_projectOut(project_out)
    {}

    template<typename GraphType, typename LinkGraphType>
    auto operator()(const GraphType& graph, const LinkGraphType& link_graph)
    {
        return run(graph, link_graph);
    }

    template<typename NodeID, typename EdgeID, typename EdgeWeight, typename LinkEdgeID, typename LinkEdgeWeight>
    std::tuple<WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>, WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>>
    run(const WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>& graph,
        const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph)
    {
        // TODO: need to check that this will work, might need to make the
        // links bidirectional for the distance oracle to work correctly
        using LinkDistance = LinkEdgeWeight;
        using Clock = std::chrono::high_resolution_clock;

        struct NodePair
        {
            NodeID u;
            NodeID v;
        };
        std::vector<bool> removable = std::vector<bool>(link_graph.num_edges(), false);
        // calculate all shortest distances between all pairs of nodes in the link graph
        WeightedTableDistOracle<NodeID, LinkEdgeID, LinkDistance> distance_oracle(link_graph);
        auto distance_func = [&](NodeID u, NodeID v) {
            return distance_oracle.get_distance(u, v);
        };

        auto block_tree_start = Clock::now();
        auto [block_tree, cycle_ids] = graph.cactus_generate_block_tree(0);
        (void)cycle_ids;
        auto block_tree_end = Clock::now();

        // root tree on node 0
        auto root_tree_start = Clock::now();
        auto [parent, depth] = rooted_tree(block_tree);
        auto root_tree_end = Clock::now();

        auto node_pairwise_dist_start = Clock::now();
        // TODO: instead of vector of vectors, could have a vector of offsets and a vector of content ala CSR
        std::vector<std::vector<NodePair>> node_pairs_by_tree_distance(graph.num_vertices());
        // We find all distances in the original graph between pairs of nodes skipping cycle nodes
        for (NodeID u{0}; u < graph.num_vertices(); u++)
        {
            for (NodeID v{u + 1}; v < graph.num_vertices(); v++)
            {
                const auto dist = tree_distance(u, v, parent, depth, graph.num_vertices());
                node_pairs_by_tree_distance[dist].emplace_back(NodePair{u, v});
            }
        }
        auto node_pairwise_dist_end = Clock::now();

        std::vector<NodeID> path = std::vector<NodeID>(graph.num_vertices());
        auto project_in_start = Clock::now();
        // Depends on the number of paths
        if (m_projectIn)
        {
            for (size_t path_length{node_pairs_by_tree_distance.size() - 1}; path_length >= 2; path_length--)
            {
                for (const NodePair& np : node_pairs_by_tree_distance[path_length])
                {
                    const NodeID u = np.u;
                    const NodeID v = np.v;
                    tree_path(u, v, parent, depth, path_length, path, graph.num_vertices());
                    const NodeID x = path[1];
                    const NodeID y = path[path_length - 1];
                    // set d(u, y) = min(d(u,v), d(u,y))
                    // and d(v, x) = min(d(u,v), d(x,v))
                    distance_oracle.set_distance(
                        u,
                        y,
                        std::min(distance_oracle.get_distance(u, v), distance_oracle.get_distance(u, y)));
                    distance_oracle.set_distance(
                        v,
                        x,
                        std::min(distance_oracle.get_distance(u, v), distance_oracle.get_distance(x, v)));
                    // Keep simmetry, bmight be unnecessary
                    distance_oracle.set_distance(y, u, distance_oracle.get_distance(u, y));
                    distance_oracle.set_distance(x, v, distance_oracle.get_distance(v, x));
                }
            }
        }
        auto project_in_end = Clock::now();

        size_t num_removed_by_project_in = 0;
        if constexpr (RecordStatsLevel > 0)
        {
            std::vector<bool> removable_after_project_in = std::vector<bool>(link_graph.num_edges(), false);
            mark_removable_links(link_graph, distance_func, removable_after_project_in);
            num_removed_by_project_in =
                std::count(removable_after_project_in.begin(), removable_after_project_in.end(), true);
        }

        auto project_out_start = Clock::now();
        if (m_projectOut)
        {
            for (size_t path_length{2}; path_length < node_pairs_by_tree_distance.size(); path_length++)
            {
                for (const NodePair& np : node_pairs_by_tree_distance[path_length])
                {
                    const NodeID u = np.u;
                    const NodeID v = np.v;
                    tree_path(u, v, parent, depth, path_length, path, graph.num_vertices());
                    LinkDistance best = distance_oracle.get_distance(u, v);
                    for (size_t i{1}; i < path_length; i++)
                    {
                        const NodeID x = path[i];
                        best = std::min(distance_oracle.get_distance(u, x) + distance_oracle.get_distance(x, v), best);
                    }
                    distance_oracle.set_distance(u, v, best);
                    // Again, symmetry might be unnecessary
                    distance_oracle.set_distance(v, u, best);
                }
            }
        }
        auto project_out_end = Clock::now();

        size_t num_removed_by_project_out = 0;
        if constexpr (RecordStatsLevel > 0)
        {
            std::vector<bool> removable_after_project_out = std::vector<bool>(link_graph.num_edges(), false);
            mark_removable_links(link_graph, distance_func, removable_after_project_out);
            num_removed_by_project_out =
                std::count(removable_after_project_out.begin(), removable_after_project_out.end(), true);
        }

        // Mark removable links
        mark_removable_links(link_graph, distance_func, removable);
        auto new_link_graph = remove_links(link_graph, removable);

        if constexpr (RecordStatsLevel > 0)
        {
            m_metrics = StageMetrics{
                {"num_removed_links", std::to_string(num_removed_by_project_out)},
                {"num_removed_by_project_in", std::to_string(num_removed_by_project_in)},
                {"num_removed_by_project_out", std::to_string(num_removed_by_project_out - num_removed_by_project_in)},
                {
                    "block_tree_construction_time",
                    HeiConnect::tools::format_duration(
                        block_tree_start,
                        block_tree_end,
                        HeiConnect::tools::TimeUnit::Seconds),
                },
                {
                    "root_tree_construction_time",
                    HeiConnect::tools::format_duration(
                        root_tree_start,
                        root_tree_end,
                        HeiConnect::tools::TimeUnit::Seconds),
                },
                {
                    "node_pairwise_distance_computation_time",
                    HeiConnect::tools::format_duration(
                        node_pairwise_dist_start,
                        node_pairwise_dist_end,
                        HeiConnect::tools::TimeUnit::Seconds),
                },
                {
                    "project_in_time",
                    HeiConnect::tools::format_duration(
                        project_in_start,
                        project_in_end,
                        HeiConnect::tools::TimeUnit::Seconds),
                },
                {
                    "project_out_time",
                    HeiConnect::tools::format_duration(
                        project_out_start,
                        project_out_end,
                        HeiConnect::tools::TimeUnit::Seconds),
                },
            };
        }

        return {graph, new_link_graph};
    }

    std::optional<StageMetrics> emit_metrics() const
    {
        if constexpr (RecordStatsLevel > 0)
        {
            return m_metrics;
        }
        else
        {
            return std::nullopt;
        }
    }

private:
    template<typename NodeID, typename EdgeID, typename WeightType, typename DistanceFn>
    void mark_removable_links(
        const WeightedCRFGraph<NodeID, EdgeID, WeightType>& link_graph,
        DistanceFn&& distance_func,
        std::vector<bool>& removable) const
    {
        for (NodeID u{0}; u < link_graph.num_vertices(); u++)
        {
            for (EdgeID e{link_graph.graph.vertices[u]}; e < link_graph.graph.vertices[u + 1]; e++)
            {
                const NodeID v = link_graph.graph.edges[e];
                const WeightType w = link_graph.weights[e];
                // TODO: we are ignoring the case where w == distance_oracle.get_distance(u, v),
                // We are doing this because you would have to check if this distance didn't take into account the
                // link (u,v) itself
                if (distance_func(u, v) < w)
                {
                    removable[e] = true;
                }
            }
        }
    }

    // Find the two nodes x and y that are adjacent to u and v respectively in the tree path from u to v.
    // P u, x, ..., y, v -> return (x, y).
    // We ignore cycle nodes in the cactus cycle graph, so that
    // P u, c1, x, ..., y, c2, v -> return (x, y) where c1 and c2 are cycle nodes.
    // TODO: Under construction
    // template<typename NodeID>
    // std::pair<NodeID, NodeID>
    // next_nodes_in_tree_path(NodeID u, NodeID v, std::vector<NodeID> parent, std::vector<size_t> depth, size_t
    // cutoff)
    // {
    //     assert false;
    //     return {u, v};
    // }

    template<typename NodeID>
    size_t tree_distance(
        NodeID u,
        NodeID v,
        const std::vector<NodeID>& parent,
        const std::vector<size_t>& depth,
        size_t cutoff) const
    {
        auto is_cycle_node = [&](NodeID node) {
            return static_cast<size_t>(node) >= cutoff;
        };

        size_t distance = 0;
        while (depth[u] > depth[v])
        {
            if (!is_cycle_node(u))
            {
                ++distance;
            }
            u = parent[u];
        }
        while (depth[v] > depth[u])
        {
            if (!is_cycle_node(v))
            {
                ++distance;
            }
            v = parent[v];
        }
        while (u != v)
        {
            if (!is_cycle_node(u))
            {
                ++distance;
            }
            if (!is_cycle_node(v))
            {
                ++distance;
            }
            u = parent[u];
            v = parent[v];
        }
        // Special case when they meet at a cycle node, then they will have added that distance twice
        if (is_cycle_node(u) && distance > 0)
            --distance;
        return distance;
    }

    // TODO: instead of distance, path, take iterators as input
    // Reconstructs the tree path from u to v using the parent and depth arrays. The result is stored in the path
    // vector.
    template<typename NodeID>
    void tree_path(
        NodeID u,
        NodeID v,
        const std::vector<NodeID>& parent,
        const std::vector<size_t>& depth,
        size_t distance,
        std::vector<NodeID>& path,
        size_t cutoff) const
    {
        auto is_cycle_node = [&](NodeID node) {
            return static_cast<size_t>(node) >= cutoff;
        };

        size_t l = 0;
        size_t r = distance;
        while (depth[u] > depth[v])
        {
            if (!is_cycle_node(u))
                path[l++] = u;
            u = parent[u];
        }
        while (depth[v] > depth[u])
        {
            if (!is_cycle_node(v))
                path[r--] = v;
            v = parent[v];
        }
        while (u != v)
        {
            if (!is_cycle_node(u))
                path[l++] = u;
            if (!is_cycle_node(v))
                path[r--] = v;
            u = parent[u];
            v = parent[v];
        }
        if (!is_cycle_node(u))
            path[l] = u;
    }

    template<typename NodeID, typename EdgeID, typename WeightType>
    std::tuple<std::vector<NodeID>, std::vector<size_t>>
    rooted_tree(const WeightedCRFGraph<NodeID, EdgeID, WeightType>& graph) const
    {
        std::vector<NodeID> parent(graph.num_vertices(), graph.num_vertices());
        std::vector<size_t> depth(graph.num_vertices(), 0);
        std::vector<bool> visited(graph.num_vertices(), false);
        std::vector<NodeID> stack;
        stack.push_back(0);
        parent[0] = 0;
        visited[0] = true;

        while (!stack.empty())
        {
            NodeID u = stack.back();
            stack.pop_back();
            for (EdgeID e{graph.graph.vertices[u]}; e < graph.graph.vertices[u + 1]; e++)
            {
                const NodeID v = graph.graph.edges[e];
                if (!visited[v])
                {
                    parent[v] = u;
                    depth[v] = depth[u] + 1;
                    visited[v] = true;
                    stack.push_back(v);
                }
            }
        }
        return {parent, depth};
    }

    bool m_projectIn;
    bool m_projectOut;
    std::optional<StageMetrics> m_metrics;
};
