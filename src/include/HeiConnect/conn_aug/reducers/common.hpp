#pragma once

#include <cstddef>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/data_structures/union_find.hpp"

using ConnAugLink = std::pair<size_t, size_t>;
using ConnAugWeightedLink = std::tuple<size_t, size_t, double>;
using ConnAugLinkRemap = std::map<ConnAugLink, ConnAugWeightedLink>;


template<typename NodeID>
ConnAugLink normalize_link(NodeID u, NodeID v)
{
    size_t a = static_cast<size_t>(u);
    size_t b = static_cast<size_t>(v);
    if (a > b)
    {
        std::swap(a, b);
    }
    return {a, b};
}


template<typename NodeID, typename LinkEdgeID, typename LinkEdgeWeight>
ConnAugLinkRemap make_identity_link_remap(const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph)
{
    ConnAugLinkRemap link_remap{};
    for (NodeID u{}; u < link_graph.num_vertices(); ++u)
    {
        for (LinkEdgeID e{link_graph.graph.vertices[u]}; e < link_graph.graph.vertices[u + 1]; ++e)
        {
            const NodeID v = link_graph.graph.edges[e];
            const ConnAugLink link = normalize_link(u, v);
            auto it = link_remap.find(link);
            if (it == link_remap.end() || static_cast<double>(link_graph.weights[e]) < std::get<2>(it->second))
            {
                link_remap[link] = {
                    static_cast<size_t>(u),
                    static_cast<size_t>(v),
                    static_cast<double>(link_graph.weights[e])};
            }
        }
    }
    return link_remap;
}


template<typename NodeID, typename LinkEdgeID, typename LinkEdgeWeight>
WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>
remove_links(const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph, const std::vector<bool>& removable)
{
    std::vector<NodeID> new_vertices(link_graph.num_vertices() + 1, 0);
    std::vector<LinkEdgeID> new_edges{};
    new_edges.reserve(link_graph.num_edges());
    std::vector<LinkEdgeWeight> new_weights{};
    new_weights.reserve(link_graph.num_edges());
    for (NodeID u{0}; u < link_graph.num_vertices(); u++)
    {
        for (LinkEdgeID e{link_graph.graph.vertices[u]}; e < link_graph.graph.vertices[u + 1]; e++)
        {
            if (!removable[e])
            {
                new_edges.push_back(link_graph.graph.edges[e]);
                new_weights.push_back(link_graph.weights[e]);
            }
        }
        new_vertices[u + 1] = new_edges.size();
    }
    return {{new_vertices, new_edges}, new_weights};
}


template<typename NodeID>
size_t
tree_distance(NodeID u, NodeID v, const std::vector<NodeID>& parent, const std::vector<size_t>& depth, size_t cutoff)
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
    size_t cutoff)
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

struct AddedLinksMergeStats
{
    size_t total_merged_nodes{0};
};

template<typename NodeID, typename Func>
void for_node_in_tree_path(
    NodeID u,
    NodeID v,
    const std::vector<NodeID>& parent,
    const std::vector<size_t>& depth,
    size_t cutoff,
    Func&& func)
{
    auto is_cycle_node = [&](NodeID node) {
        return static_cast<size_t>(node) >= cutoff;
    };

    while (depth[u] > depth[v])
    {
        if (!is_cycle_node(u))
        {
            func(u);
        }
        u = parent[u];
    }
    while (depth[v] > depth[u])
    {
        if (!is_cycle_node(v))
        {
            func(v);
        }
        v = parent[v];
    }
    while (u != v)
    {
        if (!is_cycle_node(u))
        {
            func(u);
        }
        if (!is_cycle_node(v))
        {
            func(v);
        }
        u = parent[u];
        v = parent[v];
    }
    if (!is_cycle_node(u))
    {
        func(u);
    }
}

template<typename NodeID, typename LinkEdgeID>
std::vector<NodeID> link_edge_sources(const CRFGraph<NodeID, LinkEdgeID>& link_graph)
{
    std::vector<NodeID> sources(link_graph.edges.size());
    for (NodeID u{0}; u < link_graph.num_vertices(); ++u)
    {
        for (LinkEdgeID e{link_graph.vertices[u]}; e < link_graph.vertices[u + 1]; ++e)
        {
            sources[static_cast<size_t>(e)] = u;
        }
    }
    return sources;
}

template<typename NodeID, typename LinkEdgeID, typename UnionFindType, typename LinkIDRange>
AddedLinksMergeStats add_links_to_union_find(
    size_t original_num_vertices,
    const CRFGraph<NodeID, LinkEdgeID>& link_graph,
    UnionFindType& uf,
    const LinkIDRange& link_ids,
    const std::vector<NodeID>& parent,
    const std::vector<size_t>& depth)
{
    AddedLinksMergeStats stats{};
    const auto sources = link_edge_sources(link_graph);

    for (const auto& link_id : link_ids)
    {
        const auto edge_id = static_cast<LinkEdgeID>(link_id);
        const NodeID u = sources[static_cast<size_t>(edge_id)];
        const NodeID v = link_graph.edges[edge_id];
        const size_t anchor = static_cast<size_t>(u);

        for_node_in_tree_path(u, v, parent, depth, original_num_vertices, [&](NodeID node) {
            const size_t node_id = static_cast<size_t>(node);
            if (uf.find(anchor) != uf.find(node_id))
            {
                ++stats.total_merged_nodes;
            }
            uf.unite(anchor, node_id);
        });
    }

    return stats;
}

template<typename GraphType, typename LinkGraphType, typename UnionFindType, typename LinkIDRange>
AddedLinksMergeStats add_links_to_union_find(
    const GraphType& graph,
    const LinkGraphType& link_graph,
    UnionFindType& uf,
    const LinkIDRange& link_ids)
{
    auto [block_tree, cycle_ids] = graph.cactus_generate_block_tree(0);
    (void)cycle_ids;
    auto [parent, depth] = block_tree.graph.rooted_parent_depth();
    return add_links_to_union_find(graph.num_vertices(), link_graph.graph, uf, link_ids, parent, depth);
}

template<typename NodeID, typename LinkEdgeID, typename LinkEdgeWeight, typename UnionFindType>
WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight> remap_and_eliminate_links(
    const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph,
    UnionFindType& uf,
    ConnAugLinkRemap& link_remap)
{
    ConnAugLinkRemap new_link_remap{};
    for (NodeID u{}; u < link_graph.num_vertices(); u++)
    {
        for (LinkEdgeID e{link_graph.graph.vertices[u]}; e < link_graph.graph.vertices[u + 1]; e++)
        {
            NodeID v = link_graph.graph.edges[e];
            const ConnAugLink old_link = normalize_link(u, v);
            auto weighted_link = ConnAugWeightedLink{
                static_cast<size_t>(u),
                static_cast<size_t>(v),
                static_cast<double>(link_graph.weights[e])};
            auto old_it = link_remap.find(old_link);
            if (old_it != link_remap.end())
            {
                weighted_link = old_it->second;
            }

            size_t new_u = uf.find(old_link.first);
            size_t new_v = uf.find(old_link.second);
            if (new_u != new_v)
            {
                if (new_u > new_v)
                {
                    std::swap(new_u, new_v);
                }
                ConnAugLink new_link{new_u, new_v};
                auto it = new_link_remap.find(new_link);
                if (it == new_link_remap.end() || std::get<2>(weighted_link) < std::get<2>(it->second))
                {
                    new_link_remap[new_link] = weighted_link;
                }
            }
        }
    }
    link_remap = std::move(new_link_remap);

    std::vector<NodeID> new_vertices(link_graph.num_vertices() + 1, 0);
    std::vector<LinkEdgeID> new_edges(link_remap.size());
    std::vector<LinkEdgeWeight> new_weights(link_remap.size());
    for (const auto& [new_link, weighted_link] : link_remap)
    {
        (void)weighted_link;
        ++new_vertices[new_link.first + 1];
    }
    for (size_t i{1}; i < new_vertices.size(); ++i)
    {
        new_vertices[i] += new_vertices[i - 1];
    }

    std::vector<NodeID> offsets = new_vertices;
    for (const auto& [new_link, weighted_link] : link_remap)
    {
        const auto pos = offsets[new_link.first]++;
        new_edges[pos] = static_cast<LinkEdgeID>(new_link.second);
        new_weights[pos] = static_cast<LinkEdgeWeight>(std::get<2>(weighted_link));
    }

    return {{new_vertices, new_edges}, new_weights};
}
