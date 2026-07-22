#pragma once

#include <algorithm>
#include <cstddef>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/data_structures/union_find.hpp"

using ConnAugLink = std::pair<size_t, size_t>;
struct ConnAugWeightedLink
{
    size_t original_id;
    size_t u;
    size_t v;
    double weight;
};
using ConnAugLinkRemap = std::map<ConnAugLink, ConnAugWeightedLink>;

struct ConnectivityAugmentationReductionConfig
{
    bool run_shortest_path_reduction{true};
    bool run_project_in{true};
    bool run_project_out{true};
    bool run_cycle_reduction{true};
    size_t max_rounds{0};
    bool compute_shortest_paths{true};
};


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
            if (it == link_remap.end() || static_cast<double>(link_graph.weights[e]) < it->second.weight)
            {
                link_remap[link] = {
                    static_cast<size_t>(e),
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

template<typename UnionFindType, typename CyclePositions>
size_t add_cycle_merges_to_union_find(
    size_t original_num_vertices,
    UnionFindType& uf,
    const CyclePositions& cycle_positions)
{
    std::vector<std::vector<size_t>> cycles;
    cycles.reserve(cycle_positions.size());
    for (const auto& positions : cycle_positions)
    {
        size_t cycle_size = 0;
        for (const auto position : positions)
        {
            if (position >= 0)
            {
                cycle_size = std::max(cycle_size, static_cast<size_t>(position) + 1);
            }
        }

        std::vector<size_t> cycle(cycle_size);
        for (size_t node = 0; node < positions.size(); ++node)
        {
            if (positions[node] >= 0)
            {
                cycle[static_cast<size_t>(positions[node])] = node;
            }
        }
        cycles.emplace_back(std::move(cycle));
    }

    size_t total_merged_nodes = 0;
    size_t merged_nodes = 0;
    std::vector<size_t> label_to_local(original_num_vertices, original_num_vertices);
    do
    {
        merged_nodes = 0;
        for (const auto& cycle : cycles)
        {
            std::vector<size_t> labels(cycle.size());
            std::vector<size_t> original_labels;
            original_labels.reserve(cycle.size());
            for (size_t i = 0; i < cycle.size(); ++i)
            {
                const size_t original_label = uf.find(cycle[i]);
                if (label_to_local[original_label] == original_num_vertices)
                {
                    label_to_local[original_label] = original_labels.size();
                    original_labels.emplace_back(original_label);
                }
                labels[i] = label_to_local[original_label];
            }

            std::vector<size_t> last(original_labels.size(), cycle.size());
            for (size_t i = 0; i < cycle.size(); ++i)
            {
                last[labels[i]] = i;
            }

            UnionFind cycle_uf(original_labels.size());
            std::vector<size_t> component_last = last;
            std::vector<size_t> active_index(original_labels.size(), cycle.size());
            std::vector<size_t> stack;
            stack.reserve(cycle.size());

            // Crossing components alternate on the cycle. Merge an active component with all
            // components above it to construct the smallest non-crossing partition.
            for (size_t i = 0; i < cycle.size(); ++i)
            {
                while (!stack.empty() && component_last[cycle_uf.find(stack.back())] < i)
                {
                    active_index[cycle_uf.find(stack.back())] = cycle.size();
                    stack.pop_back();
                }

                size_t label = cycle_uf.find(labels[i]);
                const size_t index = active_index[label];
                if (index == cycle.size())
                {
                    active_index[label] = stack.size();
                    stack.emplace_back(label);
                    continue;
                }
                if (index + 1 == stack.size())
                {
                    continue;
                }

                size_t new_last = component_last[label];
                for (size_t j = index + 1; j < stack.size(); ++j)
                {
                    const size_t other = cycle_uf.find(stack[j]);
                    new_last = std::max(new_last, component_last[other]);
                    active_index[other] = cycle.size();
                    cycle_uf.unite(label, other);
                    label = cycle_uf.find(label);
                }
                stack.resize(index);
                component_last[label] = new_last;
                active_index[label] = index;
                stack.emplace_back(label);
            }

            for (const size_t label : labels)
            {
                const size_t merged_label = cycle_uf.find(label);
                const size_t original_label = original_labels[label];
                const size_t original_merged_label = original_labels[merged_label];
                if (uf.find(original_label) != uf.find(original_merged_label))
                {
                    uf.unite(original_label, original_merged_label);
                    ++merged_nodes;
                }
            }
            for (const size_t original_label : original_labels)
            {
                label_to_local[original_label] = original_num_vertices;
            }
        }
        total_merged_nodes += merged_nodes;
    } while (merged_nodes > 0);

    return total_merged_nodes;
}

template<typename UnionFindType, typename CyclePositions>
size_t add_cycle_merges_to_union_find_frozen_stack(
    size_t original_num_vertices,
    UnionFindType& uf,
    const CyclePositions& cycle_positions)
{
    std::vector<size_t> base(original_num_vertices);
    for (size_t node = 0; node < original_num_vertices; ++node)
    {
        base[node] = uf.find(node);
    }

    size_t total_merged_nodes = 0;
    std::vector<size_t> label_to_local(original_num_vertices, original_num_vertices);
    for (const auto& positions : cycle_positions)
    {
        size_t cycle_size = 0;
        for (const auto position : positions)
        {
            if (position >= 0)
            {
                cycle_size = std::max(cycle_size, static_cast<size_t>(position) + 1);
            }
        }

        std::vector<size_t> labels(cycle_size);
        std::vector<size_t> original_labels;
        original_labels.reserve(cycle_size);
        for (size_t node = 0; node < positions.size(); ++node)
        {
            if (positions[node] < 0)
            {
                continue;
            }

            const size_t original_label = base[node];
            if (label_to_local[original_label] == original_num_vertices)
            {
                label_to_local[original_label] = original_labels.size();
                original_labels.emplace_back(original_label);
            }
            labels[static_cast<size_t>(positions[node])] = label_to_local[original_label];
        }

        std::vector<size_t> component_last(original_labels.size(), cycle_size);
        for (size_t i = 0; i < labels.size(); ++i)
        {
            component_last[labels[i]] = i;
        }

        UnionFind cycle_uf(original_labels.size());
        std::vector<size_t> active_index(original_labels.size(), cycle_size);
        std::vector<size_t> stack;
        stack.reserve(cycle_size);

        // Crossing components alternate on the cycle. Merge an active component with all
        // components above it to construct the smallest non-crossing partition.
        for (size_t i = 0; i < labels.size(); ++i)
        {
            while (!stack.empty() && component_last[cycle_uf.find(stack.back())] < i)
            {
                active_index[cycle_uf.find(stack.back())] = cycle_size;
                stack.pop_back();
            }

            size_t label = cycle_uf.find(labels[i]);
            const size_t index = active_index[label];
            if (index == cycle_size)
            {
                active_index[label] = stack.size();
                stack.emplace_back(label);
                continue;
            }
            if (index + 1 == stack.size())
            {
                continue;
            }

            size_t new_last = component_last[label];
            for (size_t j = index + 1; j < stack.size(); ++j)
            {
                const size_t other = cycle_uf.find(stack[j]);
                new_last = std::max(new_last, component_last[other]);
                active_index[other] = cycle_size;

                const size_t original_label = original_labels[label];
                const size_t original_other = original_labels[other];
                if (uf.find(original_label) != uf.find(original_other))
                {
                    uf.unite(original_label, original_other);
                    ++total_merged_nodes;
                }

                cycle_uf.unite(label, other);
                label = cycle_uf.find(label);
            }
            stack.resize(index);
            component_last[label] = new_last;
            active_index[label] = index;
            stack.emplace_back(label);
        }

        for (const size_t original_label : original_labels)
        {
            label_to_local[original_label] = original_num_vertices;
        }
    }

    return total_merged_nodes;
}

template<typename NodeID, typename LinkEdgeID, typename UnionFindType, typename LinkIDRange>
AddedLinksMergeStats add_links_to_union_find(
    size_t original_num_vertices,
    const CRFGraph<NodeID, LinkEdgeID>& link_graph,
    UnionFindType& uf,
    const LinkIDRange& link_ids,
    const std::vector<NodeID>& parent,
    const std::vector<size_t>& depth,
    const std::vector<std::vector<int>>& cycle_positions)
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

    stats.total_merged_nodes +=
        add_cycle_merges_to_union_find(original_num_vertices, uf, cycle_positions);

    return stats;
}

template<typename NodeID, typename LinkEdgeID, typename UnionFindType, typename LinkIDRange>
AddedLinksMergeStats add_links_to_union_find_frozen_stack(
    size_t original_num_vertices,
    const CRFGraph<NodeID, LinkEdgeID>& link_graph,
    UnionFindType& uf,
    const LinkIDRange& link_ids,
    const std::vector<NodeID>& parent,
    const std::vector<size_t>& depth,
    const std::vector<std::vector<int>>& cycle_positions)
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

    stats.total_merged_nodes +=
        add_cycle_merges_to_union_find_frozen_stack(original_num_vertices, uf, cycle_positions);

    return stats;
}

template<typename GraphType, typename LinkGraphType, typename UnionFindType, typename LinkIDRange>
AddedLinksMergeStats add_links_to_union_find(
    const GraphType& graph,
    const LinkGraphType& link_graph,
    UnionFindType& uf,
    const LinkIDRange& link_ids)
{
    auto [block_tree, cycle_positions] = graph.cactus_generate_block_tree(0);
    auto [parent, depth] = block_tree.graph.rooted_parent_depth();
    return add_links_to_union_find(
        graph.num_vertices(),
        link_graph.graph,
        uf,
        link_ids,
        parent,
        depth,
        cycle_positions);
}

template<typename GraphType, typename LinkGraphType, typename UnionFindType, typename LinkIDRange>
AddedLinksMergeStats add_links_to_union_find_frozen_stack(
    const GraphType& graph,
    const LinkGraphType& link_graph,
    UnionFindType& uf,
    const LinkIDRange& link_ids)
{
    auto [block_tree, cycle_positions] = graph.cactus_generate_block_tree(0);
    auto [parent, depth] = block_tree.graph.rooted_parent_depth();
    return add_links_to_union_find_frozen_stack(
        graph.num_vertices(),
        link_graph.graph,
        uf,
        link_ids,
        parent,
        depth,
        cycle_positions);
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
                static_cast<size_t>(e),
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
                if (it == new_link_remap.end() || weighted_link.weight < it->second.weight)
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
        new_weights[pos] = static_cast<LinkEdgeWeight>(weighted_link.weight);
    }

    return {{new_vertices, new_edges}, new_weights};
}

template<typename NodeID, typename EdgeID, typename EdgeWeight, typename LinkEdgeID, typename LinkEdgeWeight>
auto materialize_contractions(
    const WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>& graph,
    const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph,
    UnionFind& uf,
    ConnAugLinkRemap& link_remap)
{
    std::map<size_t, size_t> contracted_ids;
    std::vector<size_t> node_remap(graph.num_vertices());
    for (size_t u = 0; u < graph.num_vertices(); ++u)
    {
        const size_t representative = uf.find(u);
        auto [it, inserted] = contracted_ids.try_emplace(representative, contracted_ids.size());
        (void)inserted;
        node_remap[u] = it->second;
    }

    std::map<ConnAugLink, EdgeWeight> contracted_edges;
    for (NodeID u{}; u < graph.num_vertices(); ++u)
    {
        for (EdgeID e = graph.graph.vertices[u]; e < graph.graph.vertices[u + 1]; ++e)
        {
            const NodeID v = graph.graph.edges[e];
            if (u > v)
            {
                continue;
            }
            const ConnAugLink edge = normalize_link(node_remap[u], node_remap[v]);
            if (edge.first != edge.second)
            {
                contracted_edges[edge] += graph.weights[e];
            }
        }
    }

    std::vector<std::tuple<NodeID, NodeID, EdgeWeight>> graph_edges;
    graph_edges.reserve(contracted_edges.size() * 2);
    for (const auto& [edge, weight] : contracted_edges)
    {
        graph_edges.emplace_back(edge.first, edge.second, weight);
        graph_edges.emplace_back(edge.second, edge.first, weight);
    }
    auto contracted_graph = WeightedCRFGraph<NodeID, EdgeID, EdgeWeight>::vec_links_to_csr(
        graph_edges,
        contracted_ids.size());

    ConnAugLinkRemap contracted_link_remap;
    for (NodeID u{}; u < link_graph.num_vertices(); ++u)
    {
        for (LinkEdgeID e = link_graph.graph.vertices[u]; e < link_graph.graph.vertices[u + 1]; ++e)
        {
            const NodeID v = link_graph.graph.edges[e];
            const auto old_it = link_remap.find(normalize_link(u, v));
            if (old_it == link_remap.end())
            {
                continue;
            }
            const ConnAugLink link = normalize_link(node_remap[u], node_remap[v]);
            if (link.first == link.second)
            {
                continue;
            }
            auto it = contracted_link_remap.find(link);
            if (it == contracted_link_remap.end() || old_it->second.weight < it->second.weight)
            {
                contracted_link_remap[link] = old_it->second;
            }
        }
    }

    std::vector<std::tuple<NodeID, NodeID, LinkEdgeWeight>> links;
    links.reserve(contracted_link_remap.size());
    for (const auto& [link, origin] : contracted_link_remap)
    {
        links.emplace_back(link.first, link.second, static_cast<LinkEdgeWeight>(origin.weight));
    }
    auto contracted_link_graph = WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>::vec_links_to_csr(
        links,
        contracted_ids.size());
    link_remap = std::move(contracted_link_remap);

    return std::tuple{
        std::move(contracted_graph),
        std::move(contracted_link_graph),
        std::move(node_remap)};
}
