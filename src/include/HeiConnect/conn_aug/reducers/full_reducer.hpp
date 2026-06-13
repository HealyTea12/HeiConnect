#pragma once

#include "HeiConnect/data_structures/distance_oracle.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/conn_aug/reducers/common.hpp"

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
class FullReducer
{
public:
    FullReducer(bool project_in, bool project_out) : m_projectIn(project_in), m_projectOut(project_out) {};
    template<typename GraphType, typename LinkGraphType>
    auto operator()(const GraphType& graph, const LinkGraphType& link_graph)
    {
        return run(graph, link_graph);
    }
    template<typename GraphType, typename LinkGraphType>
    auto run(const GraphType& graph, const LinkGraphType& link_graph)
    {
        // TODO: need to check that this will work, might need to make the
        // links bidirectional for the distance oracle to work correctly
        using NodeID = typename GraphType::NodeID;
        using EdgeID = typename GraphType::EdgeID;
        using Distance = typename GraphType::EdgeWeight;
        std::vector<bool> removable = std::vector<bool>(link_graph.num_edges(), false);
        // calculate all shortest distances between all pairs of nodes in the link graph
        WeightedTableDistOracle<NodeID, EdgeID, Distance> distance_oracle(link_graph);

        // root tree on node 0
        auto [parent, depth] = rooted_tree(graph);
        if (m_projectIn)
        {
            // we need to do this for every link, so it doesn't have to be bidirectional
            // for (NodeID u{0}; u < link_graph.num_vertices(); u++)
            // {
            //     for (auto [v, w] : link_graph.get_neighbors(u))
            //     {
            //         project_in(u, v, distance_oracle, parent, depth);
            //     }
            // }
            // Find all leaves in the original graph
            std::vector<NodeID> leaves;
            leaves.reserve(graph.num_vertices());
            graph.forEachVertex([&](NodeID u) {
                if (graph.degree(u) == 1)
                {
                    leaves.emplace_back(u);
                }
            });
            // Propagate distances inwards from the leaves towards the root
            std::vector<bool> visited(graph.num_vertices(), false);
            std::vector<NodeID> new_leaves{};
            new_leaves.reserve(graph.num_vertices());
            while (!leaves.empty())
            {
                new_leaves.clear();
                for (size_t leaf_1_idx{}; leaf_1_idx < leaves.size(); leaf_1_idx++)
                {
                    NodeID leaf_1 = leaves[leaf_1_idx];
                    for (size_t leaf_2_idx{leaf_1_idx + 1}; leaf_2_idx < leaves.size(); leaf_2_idx++)
                    {
                        NodeID leaf_2 = leaves[leaf_2_idx];
                        // add the parents to the next batch of leaves to process
                        NodeID parent_1 = parent[leaf_1];
                        NodeID parent_2 = parent[leaf_2];
                        if (!visited[parent_1])
                        {
                            new_leaves.emplace_back(parent_1);
                            visited[parent_1] = true;
                        }
                        if (!visited[parent_2])
                        {
                            new_leaves.emplace_back(parent_2);
                            visited[parent_2] = true;
                        }
                        // project distance d(leaf_1, leaf_2) to (parent_1, leaf_2), (leaf_1, parent_2) and (parent_1,
                        // parent_2)
                        auto d_l1_l2 = distance_oracle.get_distance(leaf_1, leaf_2);
                        auto d_p1_l2 = distance_oracle.get_distance(parent_1, leaf_2);
                        auto d_l1_p2 = distance_oracle.get_distance(leaf_1, parent_2);
                        auto d_p1_p2 = distance_oracle.get_distance(parent_1, parent_2);
                        distance_oracle.set_distance(parent_1, leaf_2, std::min(d_p1_l2, d_l1_l2));
                        distance_oracle.set_distance(leaf_2, parent_1, std::min(d_p1_l2, d_l1_l2));
                        distance_oracle.set_distance(leaf_1, parent_2, std::min(d_l1_p2, d_l1_l2));
                        distance_oracle.set_distance(parent_2, leaf_1, std::min(d_l1_p2, d_l1_l2));
                        if (parent_1 != parent_2)
                        {
                            distance_oracle.set_distance(
                                parent_1,
                                parent_2,
                                std::min({d_p1_p2, d_p1_l2, d_l1_p2, d_l1_l2}));
                            distance_oracle.set_distance(
                                parent_2,
                                parent_1,
                                std::min({d_p1_p2, d_p1_l2, d_l1_p2, d_l1_l2}));
                        }
                    }
                }
                std::swap(leaves, new_leaves);
            }
        }
        if (m_projectOut)
        {
            struct NodePair
            {
                NodeID u;
                NodeID v;
            };
            std::vector<std::vector<NodePair>> node_pairs_by_tree_distance(graph.num_vertices());
            // Original graph (og)
            TableDistOracle og_distance_oracle(graph);
            for (NodeID u{0}; u < graph.num_vertices(); u++)
            {
                for (NodeID v{u + 1}; v < graph.num_vertices(); v++)
                {
                    size_t tree_distance = og_distance_oracle.get_distance(u, v);
                    node_pairs_by_tree_distance[tree_distance].emplace_back(NodePair{u, v});
                }
            }
            std::vector<NodeID> path = std::vector<NodeID>(node_pairs_by_tree_distance.size() + 1);
            for (size_t path_length{2}; path_length < node_pairs_by_tree_distance.size(); path_length++)
            {
                for (const NodePair& np : node_pairs_by_tree_distance[path_length])
                {
                    NodeID u = np.u;
                    NodeID v = np.v;
                    tree_path(u, v, parent, depth, path, path_length);
                    for (size_t i{1}; i < path_length; i++)
                    {
                        NodeID x = path[i];
                        Distance new_distance = std::min(
                            distance_oracle.get_distance(u, x) + distance_oracle.get_distance(x, v),
                            distance_oracle.get_distance(u, v));
                    }
                }
            }
        }

        for (NodeID u{0}; u < link_graph.num_vertices(); u++)
        {
            for (const auto& [v, w] : link_graph.get_neighbors(u))
            {
                if (distance_oracle.get_distance(u, v) < w)
                {
                    removable[link_graph.edge_id(u, v)] = true;
                }
            }
        }
        LinkGraphType new_link_graph = remove_links(link_graph, removable);
        return {graph, new_link_graph};
    }

private:
    // Reconstructs the tree path from u to v using the parent and depth arrays. The result is stored in the path
    // vector.
    template<typename NodeID>
    void tree_path(
        NodeID u,
        NodeID v,
        std::vector<NodeID> parent,
        std::vector<size_t> depth,
        size_t distance,
        std::vector<NodeID>& path)
    {
        size_t l = 0;
        size_t r = distance;
        while (u != v)
        {
            if (depth[u] < depth[v])
            {
                path[r--] = v;
                v = parent[v];
            }
            else
            {
                path[l++] = u;
                u = parent[u];
            }
        }
        path[l] = u;
    }


    // The two links closest to u<->v on the path from u to v get assigned the minimum of
    // the two distances to u and v and themselves
    void project_in(
        size_t u,
        size_t v,
        WeightedTableDistOracle<size_t>& distance_oracle,
        const std::vector<size_t>& parent,
        const std::vector<size_t>& depth)
    {
        size_t v_next;
        size_t u_next;
    }
    template<typename GraphType>
    std::tuple<std::vector<typename GraphType::NodeID>, std::vector<size_t>> rooted_tree(const GraphType& graph)
    {
        using NodeID = typename GraphType::NodeID;
        std::vector<NodeID> parent(graph.num_vertices(), graph.num_vertices());
        std::vector<size_t> depth(graph.num_vertices(), 0);
        std::vector<bool> visited(graph.num_vertices(), false);
        std::vector<NodeID> stack;
        stack.push_back(0);
        parent[0] = 0;
        while (!stack.empty())
        {
            NodeID u = stack.back();
            stack.pop_back();
            visited[u] = true;
            for (const auto& [v, w] : graph.get_neighbors(u))
            {
                if (!visited[v])
                {
                    parent[v] = u;
                    depth[v] = depth[u] + 1;
                    stack.push_back(v);
                }
            }
        }
        return {parent, depth};
    }

private:
    bool m_projectIn;
    bool m_projectOut;
};