#pragma once

#include <algorithm>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "HeiConnect/conn_aug/reducers/common.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/data_structures/union_find.hpp"
#include "HeiConnect/pipeline/common.hpp"


/*TODO:
The idea of this reducer is that for every edge $e_dominator$ in the cycle-path between the endpoints of a link,
if whenever $e_dominator is covered by a link, an edge $e_dominee$ is also covered by the link,
then $e_dominee$ can be removed (i.e contracted) from the graph.

To do so, we need to compute two data structures.
One that gives for every link, the (sorted) ids of the edges that are covered by it and
one that gives for every edge, the ids of the links that cover it.

Then:
for each edge:
    l_small <- link of smallest degree
    candidates <- link_to_edges(l_small)
    for each link in edge_to_links(edge):
        candidates <- candidates intersect link_to_edges(link)
*/

class FullMinCutDomReducer
{
public:
    static std::string name()
    {
        return "Element Domination Reducer";
    }

    using LinkRemap = ConnAugLinkRemap;

    FullMinCutDomReducer() = default;

    template<typename NodeID, typename EdgeID, typename EdgeWeight, typename LinkEdgeID, typename LinkWeight>
    void operator()(
        const WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>& graph,
        const WeightedCRFGraph<NodeID, LinkEdgeID, LinkWeight>& link_graph,
        UnionFind& uf)
    {
        auto [dominated_edges, edge_endpoints] = collect_dominated_edges(graph, link_graph);
        for (size_t dominated_edge : dominated_edges)
        {
            const auto it = edge_endpoints.find(dominated_edge);
            if (it == edge_endpoints.end())
            {
                continue;
            }
            uf.unite(static_cast<size_t>(it->second.first), static_cast<size_t>(it->second.second));
        }
    }

    template<typename GraphType, typename LinkGraphType>
    auto operator()(const GraphType& graph, const LinkGraphType& link_graph, LinkRemap& link_remap, UnionFind& uf)
    {
        return run(graph, link_graph, link_remap, uf);
    }

    template<typename NodeID, typename EdgeID, typename EdgeWeight, typename LinkEdgeID, typename LinkWeight>
    std::tuple<
        WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>,
        WeightedCRFGraph<NodeID, LinkEdgeID, LinkWeight>,
        LinkRemap,
        UnionFind>
    run(const WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>& graph,
        const WeightedCRFGraph<NodeID, LinkEdgeID, LinkWeight>& link_graph,
        LinkRemap& link_remap,
        UnionFind& uf)
    {
        auto [dominated_edges, edge_endpoints] = collect_dominated_edges(graph, link_graph);
        size_t num_contracted_pairs = 0;
        for (size_t dominated_edge : dominated_edges)
        {
            const auto it = edge_endpoints.find(dominated_edge);
            if (it == edge_endpoints.end())
            {
                continue;
            }
            const size_t u = static_cast<size_t>(it->second.first);
            const size_t v = static_cast<size_t>(it->second.second);
            if (uf.find(u) != uf.find(v))
            {
                ++num_contracted_pairs;
            }
            uf.unite(u, v);
        }

        auto new_link_graph = remap_and_eliminate_links(link_graph, uf, link_remap);
        m_metrics = StageMetrics{
            {"num_dominated_pairs", std::to_string(dominated_edges.size())},
            {"num_contracted_pairs", std::to_string(num_contracted_pairs)},
            {"num_removed_links", std::to_string(link_graph.num_edges() - new_link_graph.num_edges())},
        };

        return {graph, new_link_graph, link_remap, uf};
    }

    std::optional<StageMetrics> emit_metrics() const
    {
        return m_metrics;
    }

private:
    template<typename NodeID, typename EdgeID, typename EdgeWeight, typename LinkEdgeID, typename LinkWeight>
    std::pair<std::vector<size_t>, std::unordered_map<size_t, std::pair<NodeID, NodeID>>> collect_dominated_edges(
        const WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>& graph,
        const WeightedCRFGraph<NodeID, LinkEdgeID, LinkWeight>& link_graph)
    {
        std::vector<std::tuple<NodeID, NodeID, LinkWeight>> links = link_graph.csr_to_vec_links();
        auto [block_tree, cycle_positions] = graph.cactus_generate_block_tree(0);
        auto [parent, depth] = block_tree.graph.rooted_parent_depth();

        std::unordered_map<size_t, std::pair<NodeID, NodeID>> edge_endpoints{};
        std::vector<size_t> cycle_size(cycle_positions.size());
        for (size_t cid = 0; cid < cycle_positions.size(); ++cid)
        {
            for (const auto pos : cycle_positions[cid])
            {
                if (pos >= 0)
                {
                    cycle_size[cid] = static_cast<size_t>(pos) + 1;
                }
            }
        }
        auto link_to_edges =
            construct_links_to_edges(graph, links, parent, depth, cycle_positions, cycle_size, edge_endpoints);
        auto edge_to_links = construct_edges_to_links(link_to_edges);
        return {find_dominated_edges(link_to_edges, edge_to_links), std::move(edge_endpoints)};
    }

    template<typename NodeID>
    size_t edge_id(NodeID u, NodeID v, size_t num_vertices) const
    {
        if (u > v)
        {
            std::swap(u, v);
        }
        return static_cast<size_t>(u) * num_vertices + static_cast<size_t>(v);
    }

    template<typename NodeID, typename EdgeID, typename EdgeWeight, typename LinkWeight>
    std::vector<std::vector<size_t>> construct_links_to_edges(
        const WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>& graph,
        const std::vector<std::tuple<NodeID, NodeID, LinkWeight>>& links,
        const std::vector<NodeID>& parent,
        const std::vector<size_t>& depth,
        const std::vector<std::vector<int>>& cycle_positions,
        const std::vector<size_t>& cycle_size,
        std::unordered_map<size_t, std::pair<NodeID, NodeID>>& edge_endpoints) const
    {
        std::vector<std::vector<size_t>> link_to_edges(links.size());
        for (size_t lid = 0; lid < links.size(); ++lid)
        {
            const auto [u, v, w] = links[lid];
            (void)w;
            const auto original_num_vertices = graph.num_vertices();
            auto is_cycle_node = [&](NodeID node) {
                return static_cast<size_t>(node) >= original_num_vertices;
            };

            std::vector<NodeID> path_u_to_lca;
            std::vector<NodeID> path_lca_to_v;
            NodeID uu = u;
            NodeID vv = v;
            while (depth[uu] > depth[vv])
            {
                path_u_to_lca.emplace_back(uu);
                uu = parent[uu];
            }
            while (depth[vv] > depth[uu])
            {
                path_lca_to_v.emplace_back(vv);
                vv = parent[vv];
            }
            while (uu != vv)
            {
                path_u_to_lca.emplace_back(uu);
                path_lca_to_v.emplace_back(vv);
                uu = parent[uu];
                vv = parent[vv];
            }
            path_u_to_lca.emplace_back(uu);
            std::reverse(path_lca_to_v.begin(), path_lca_to_v.end());
            path_u_to_lca.insert(path_u_to_lca.end(), path_lca_to_v.begin(), path_lca_to_v.end());

            auto& path = path_u_to_lca;

            auto& covered_edges = link_to_edges[lid];
            covered_edges.reserve(path.size() == 0 ? 0 : path.size() - 1);
            for (size_t i = 0; i + 1 < path.size(); ++i)
            {
                NodeID a = path[i];
                NodeID b = path[i + 1];
                if (is_cycle_node(a) || is_cycle_node(b))
                {
                    continue;
                }
                if (a > b)
                {
                    std::swap(a, b);
                }
                const size_t id = edge_id(a, b, original_num_vertices);
                covered_edges.emplace_back(id);
                edge_endpoints.emplace(id, std::make_pair(a, b));
            }

            for (size_t i = 0; i + 2 < path.size(); ++i)
            {
                if (!is_cycle_node(path[i]) && is_cycle_node(path[i + 1]) && !is_cycle_node(path[i + 2]))
                {
                    NodeID a = path[i];
                    NodeID b = path[i + 2];
                    const auto cycle_id = static_cast<size_t>(path[i + 1]) - original_num_vertices;
                    if (cycle_id >= cycle_positions.size())
                    {
                        continue;
                    }
                    if (cycle_size[cycle_id] < 4 || cycle_size[cycle_id] % 2 != 0)
                    {
                        continue;
                    }
                    const int cycle_pos_a = cycle_positions[cycle_id][static_cast<size_t>(a)];
                    const int cycle_pos_b = cycle_positions[cycle_id][static_cast<size_t>(b)];
                    if (cycle_pos_a < 0 || cycle_pos_b < 0)
                    {
                        continue;
                    }
                    const auto pos_a = static_cast<size_t>(cycle_pos_a);
                    const auto pos_b = static_cast<size_t>(cycle_pos_b);
                    const auto gap = (pos_a > pos_b) ? (pos_a - pos_b) : (pos_b - pos_a);
                    const auto cycle_gap = std::min(gap, cycle_size[cycle_id] - gap);
                    if (cycle_gap != cycle_size[cycle_id] / 2)
                    {
                        continue;
                    }
                    if (a < b)
                    {
                        const size_t id = edge_id(a, b, original_num_vertices);
                        covered_edges.emplace_back(id);
                        edge_endpoints.emplace(id, std::make_pair(a, b));
                    }
                    else
                    {
                        const size_t id = edge_id(b, a, original_num_vertices);
                        covered_edges.emplace_back(id);
                        edge_endpoints.emplace(id, std::make_pair(b, a));
                    }
                }
            }

            std::sort(covered_edges.begin(), covered_edges.end());
            covered_edges.erase(std::unique(covered_edges.begin(), covered_edges.end()), covered_edges.end());
        }
        return link_to_edges;
    }

    std::unordered_map<size_t, std::vector<size_t>>
    construct_edges_to_links(const std::vector<std::vector<size_t>>& link_to_edges) const
    {
        std::unordered_map<size_t, std::vector<size_t>> edge_to_links{};
        for (size_t lid = 0; lid < link_to_edges.size(); ++lid)
        {
            for (size_t covered_edge : link_to_edges[lid])
            {
                edge_to_links[covered_edge].emplace_back(lid);
            }
        }
        return edge_to_links;
    }

    std::vector<size_t> find_dominated_edges(
        const std::vector<std::vector<size_t>>& link_to_edges,
        const std::unordered_map<size_t, std::vector<size_t>>& edge_to_links) const
    {
        std::vector<size_t> dominated_edges{};
        for (const auto& [dominator_edge, covering_links] : edge_to_links)
        {
            if (covering_links.empty())
            {
                continue;
            }

            size_t smallest_link = covering_links.front();
            for (size_t lid : covering_links)
            {
                if (link_to_edges[lid].size() < link_to_edges[smallest_link].size())
                {
                    smallest_link = lid;
                }
            }

            std::vector<size_t> candidates = link_to_edges[smallest_link];
            for (size_t lid : covering_links)
            {
                if (lid == smallest_link)
                {
                    continue;
                }

                std::vector<size_t> intersection{};
                intersection.reserve(std::min(candidates.size(), link_to_edges[lid].size()));
                std::set_intersection(
                    candidates.begin(),
                    candidates.end(),
                    link_to_edges[lid].begin(),
                    link_to_edges[lid].end(),
                    std::back_inserter(intersection));
                candidates = std::move(intersection);
                if (candidates.empty())
                {
                    break;
                }
            }

            for (size_t candidate_edge : candidates)
            {
                if (candidate_edge != dominator_edge)
                {
                    dominated_edges.emplace_back(candidate_edge);
                }
            }
        }

        std::sort(dominated_edges.begin(), dominated_edges.end());
        dominated_edges.erase(std::unique(dominated_edges.begin(), dominated_edges.end()), dominated_edges.end());
        return dominated_edges;
    }
    std::optional<StageMetrics> m_metrics;
};
