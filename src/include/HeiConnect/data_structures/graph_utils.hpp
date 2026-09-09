#pragma once
#include <algorithm>
#include <numeric>
#include <random>
#include <stdexcept>

#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/graph.hpp"

inline WeightedCRFGraph<> create_cycle_graph(size_t n_nodes)
{
    auto vertices = std::vector<size_t>(n_nodes + 1, 0);
    std::vector<size_t> edges{};
    std::vector<double> weights{};
    for (size_t u = 0; u < n_nodes; u++)
    {
        edges.push_back((u + 1) % n_nodes);
        weights.push_back(1.0);
        vertices[u + 1] = edges.size();
    }
    return WeightedCRFGraph<>{CRFGraph<>{vertices, edges}, weights};
}

inline WeightedCRFGraph<> create_cycle_graph_undirected(size_t n_nodes)
{
    auto vertices = std::vector<size_t>(n_nodes + 1, 0);
    std::vector<size_t> edges{};
    std::vector<double> weights{};
    for (size_t u = 0; u < n_nodes; u++)
    {
        edges.push_back((u + 1) % n_nodes);
        weights.push_back(1.0);
        edges.push_back((u + n_nodes - 1) % n_nodes);
        weights.push_back(1.0);
        vertices[u + 1] = edges.size();
    }
    return WeightedCRFGraph<>{CRFGraph<>{vertices, edges}, weights};
}

// Creates a star graph with one center node (node 0) connected to n_leaves leaf nodes.
template <class node_T = uint64_t, class edge_T = uint64_t, class weight_T = double>
inline WeightedCRFGraph<node_T, edge_T, weight_T> create_star_graph(node_T n_leaves)
{
    auto vertices = std::vector<edge_T>(n_leaves + 2, static_cast<edge_T>(0));
    auto edges = std::vector<node_T>(2 * n_leaves, static_cast<node_T>(0));
    auto weights = std::vector<weight_T>(2 * n_leaves, static_cast<weight_T>(1.0));
    vertices[1] = n_leaves;
    for (node_T i = static_cast<node_T>(1); i <= n_leaves; i++)
    {
        vertices[i + 1] = vertices[i] + 1;
        edges[i - 1] = i; // center to leaf
        edges[n_leaves + i - 1] = 0;
    }
    return WeightedCRFGraph<node_T, edge_T, weight_T>{{vertices, edges}, weights};
}

inline WeightedCRFGraph<> create_random_tree(size_t n_nodes, unsigned int seed = 42)
{
    // tree has n - 1 edges
    auto vertices = std::vector<size_t>(n_nodes + 1, 0);
    auto edges = std::vector<size_t>(2 * (n_nodes - 1), 0);
    auto weights = std::vector<double>(2 * (n_nodes - 1), 1.0);
    auto random_engine = std::mt19937(seed);
    auto temp_edges = std::vector<std::vector<size_t>>(n_nodes, std::vector<size_t>{});
    for (size_t u = 1; u < n_nodes; u++)
    {
        // choose one of the already connected nodes
        std::uniform_int_distribution<size_t> dist{0, u - 1};
        auto v = dist(random_engine);
        temp_edges[u].emplace_back(v);
        temp_edges[v].emplace_back(u);
    }
    for (size_t u = 0; u < n_nodes; u++)
    {
        auto starting_point = vertices[u]; // we start placing us neighbours here
        vertices[u + 1] = starting_point + temp_edges[u].size();
        for (size_t i = 0; i < temp_edges[u].size(); i++)
        {
            auto v = temp_edges[u][i];
            edges[starting_point + i] = v;
        }
    }
    return WeightedCRFGraph<>{{vertices, edges}, weights};
}

/**
 * Creates a connected random cactus in weighted CSR form.
 *
 * The graph starts with one node. It then adds bridge blocks and fixed-length
 * cycle blocks in random order. Every block is attached to a uniformly chosen
 * node that is already part of the graph. Finally, the node labels are randomly
 * permuted so that they do not reveal the order in which the nodes were created.
 * Cycle edges receive weight 1.0 and bridge edges receive weight 2.0, so both a
 * cut through two cycle edges and a cut through one bridge have value 2.0.
 *
 * @param node_count Total number of nodes in the generated cactus. This must be
 *        at least one and includes the initial node.
 * @param cycle_count Number of cycle blocks to add. Setting this to zero
 *        generates a random tree.
 * @param cycle_length Number of nodes in every cycle block. This must be at
 *        least three. Each cycle uses one existing attachment node and creates
 *        cycle_length - 1 new nodes.
 * @param seed Seed used for block shuffling, attachment-node selection, and the
 *        final node-label permutation. The same parameters and seed produce the
 *        same graph.
 * @return An undirected WeightedCRFGraph whose adjacency arrays store both
 *         directions of every edge.
 * @throws std::invalid_argument If node_count or cycle_length is invalid, or if
 *         cycle_count * (cycle_length - 1) exceeds node_count - 1.
 */
inline WeightedCRFGraph<> create_random_cactus(
    size_t node_count,
    size_t cycle_count,
    size_t cycle_length,
    unsigned int seed = 42)
{
    if (node_count == 0)
        throw std::invalid_argument("a cactus must contain at least one node");
    if (cycle_length < 3)
        throw std::invalid_argument("cactus cycles must contain at least three nodes");
    if (cycle_count > 0 && cycle_length - 1 > (node_count - 1) / cycle_count)
        throw std::invalid_argument("the requested cycles require more than the available nodes");

    enum class BlockType
    {
        bridge,
        cycle
    };

    const size_t cycle_node_count = cycle_count * (cycle_length - 1);
    const size_t bridge_count = node_count - 1 - cycle_node_count;

    auto random_engine = std::mt19937(seed);
    auto blocks = std::vector<BlockType>(bridge_count, BlockType::bridge);
    blocks.insert(blocks.end(), cycle_count, BlockType::cycle);
    std::shuffle(blocks.begin(), blocks.end(), random_engine);

    auto adjacency_lists = std::vector<std::vector<std::pair<size_t, double>>>(node_count);
    size_t next_node = 1;
    for (const auto block : blocks)
    {
        std::uniform_int_distribution<size_t> attachment_distribution(0, next_node - 1);
        const size_t attachment_node = attachment_distribution(random_engine);

        if (block == BlockType::bridge)
        {
            const size_t new_node = next_node++;
            adjacency_lists[attachment_node].emplace_back(new_node, 2.0);
            adjacency_lists[new_node].emplace_back(attachment_node, 2.0);
            continue;
        }

        size_t previous_node = attachment_node;
        for (size_t cycle_node = 0; cycle_node < cycle_length - 1; ++cycle_node)
        {
            const size_t new_node = next_node++;
            adjacency_lists[previous_node].emplace_back(new_node, 1.0);
            adjacency_lists[new_node].emplace_back(previous_node, 1.0);
            previous_node = new_node;
        }
        adjacency_lists[previous_node].emplace_back(attachment_node, 1.0);
        adjacency_lists[attachment_node].emplace_back(previous_node, 1.0);
    }

    auto shuffled_labels = std::vector<size_t>(node_count);
    std::iota(shuffled_labels.begin(), shuffled_labels.end(), 0);
    std::shuffle(shuffled_labels.begin(), shuffled_labels.end(), random_engine);

    auto relabeled_adjacency_lists = std::vector<std::vector<std::pair<size_t, double>>>(node_count);
    for (size_t node = 0; node < node_count; ++node)
    {
        for (const auto &[neighbor, weight] : adjacency_lists[node])
        {
            relabeled_adjacency_lists[shuffled_labels[node]].emplace_back(
                shuffled_labels[neighbor], weight);
        }
    }

    auto vertices = std::vector<size_t>(node_count + 1, 0);
    auto edges = std::vector<size_t>();
    auto weights = std::vector<double>();
    edges.reserve(2 * (node_count - 1 + cycle_count));
    weights.reserve(edges.capacity());
    for (size_t node = 0; node < node_count; ++node)
    {
        for (const auto &[neighbor, weight] : relabeled_adjacency_lists[node])
        {
            edges.push_back(neighbor);
            weights.push_back(weight);
        }
        vertices[node + 1] = edges.size();
    }

    return WeightedCRFGraph<>{{vertices, edges}, weights};
}

// Creates exactly cycle_count blocks on exactly node_count vertices, counting
// bridges as cycles of length two. Start each block with one new vertex and one
// attachment. Each surplus vertex uniformly chooses one of the blocks to extend.
// Blocks that remain of size two become bridges. Shuffle the blocks and attach each at a
// uniformly chosen existing vertex. This is not a uniform sample of all cacti.
// Labels are randomly permuted; cycle edges have weight 1 and bridges weight 2.
// Setting cycle_count to node_count - 1 produces a tree.
// Requires 1 <= cycle_count <= node_count - 1, or node_count = 1 and cycle_count = 0.
// Construction takes O(node_count) time and space and is reproducible by seed.
inline WeightedCRFGraph<> create_random_cactus_with_cycle_count(
    size_t node_count,
    size_t cycle_count,
    unsigned int seed = 42)
{
    if (node_count == 0)
        throw std::invalid_argument("a cactus must contain at least one node");
    if (cycle_count > node_count - 1)
        throw std::invalid_argument("the requested cycles require more than the available nodes");
    if (cycle_count == 0 && node_count > 1)
        throw std::invalid_argument("a cactus with more than one node requires at least one cycle or bridge");

    auto random_engine = std::mt19937(seed);
    auto block_sizes = std::vector<size_t>(cycle_count, 2);
    const size_t surplus = node_count - 1 - cycle_count;
    if (cycle_count > 0)
    {
        auto allocation = std::uniform_int_distribution<size_t>(0, cycle_count - 1);
        for (size_t node = 0; node < surplus; ++node)
            ++block_sizes[allocation(random_engine)];
    }
    std::shuffle(block_sizes.begin(), block_sizes.end(), random_engine);

    auto labels = std::vector<size_t>(node_count);
    std::iota(labels.begin(), labels.end(), 0);
    std::shuffle(labels.begin(), labels.end(), random_engine);
    auto adjacency_lists = std::vector<std::vector<std::pair<size_t, double>>>(node_count);
    auto add_edge = [&](size_t u, size_t v, double weight)
    {
        adjacency_lists[labels[u]].emplace_back(labels[v], weight);
        adjacency_lists[labels[v]].emplace_back(labels[u], weight);
    };

    size_t next_node = 1;
    for (const size_t block_size : block_sizes)
    {
        auto attachment_distribution = std::uniform_int_distribution<size_t>(0, next_node - 1);
        const size_t attachment = attachment_distribution(random_engine);
        if (block_size == 2)
        {
            add_edge(attachment, next_node++, 2.0);
            continue;
        }
        size_t previous = attachment;
        for (size_t node = 1; node < block_size; ++node)
        {
            const size_t new_node = next_node++;
            add_edge(previous, new_node, 1.0);
            previous = new_node;
        }
        add_edge(previous, attachment, 1.0);
    }

    auto vertices = std::vector<size_t>(node_count + 1, 0);
    auto edges = std::vector<size_t>();
    auto weights = std::vector<double>();
    edges.reserve(2 * (node_count - 1 + cycle_count));
    weights.reserve(edges.capacity());
    for (size_t node = 0; node < node_count; ++node)
    {
        for (const auto &[neighbor, weight] : adjacency_lists[node])
        {
            edges.push_back(neighbor);
            weights.push_back(weight);
        }
        vertices[node + 1] = edges.size();
    }
    return WeightedCRFGraph<>{{vertices, edges}, weights};
}

// Attaches cycles whose sizes are sampled uniformly from the inclusive bounds.
// A cycle uses one existing vertex, so it introduces cycle_size - 1 vertices.
// The last addition is truncated to fit; a single remaining vertex gets a bridge.
// A sampled cycle size of two also adds a bridge with one new vertex.
// Cycle edges have weight 1 and bridges have weight 2, as above.
inline WeightedCRFGraph<> create_random_cactus_with_cycle_sizes(
    size_t node_count,
    size_t min_cycle_size,
    size_t max_cycle_size,
    unsigned int seed = 42)
{
    if (node_count == 0)
        throw std::invalid_argument("a cactus must contain at least one node");
    if (min_cycle_size < 2 || max_cycle_size < min_cycle_size)
        throw std::invalid_argument("cycle sizes must satisfy 2 <= min_cycle_size <= max_cycle_size");

    auto random_engine = std::mt19937(seed);
    auto cycle_size_distribution = std::uniform_int_distribution<size_t>(min_cycle_size, max_cycle_size);
    auto adjacency_lists = std::vector<std::vector<std::pair<size_t, double>>>(node_count);
    auto add_edge = [&](size_t u, size_t v, double weight)
    {
        adjacency_lists[u].emplace_back(v, weight);
        adjacency_lists[v].emplace_back(u, weight);
    };

    size_t next_node = 1;
    while (next_node < node_count)
    {
        auto attachment_distribution = std::uniform_int_distribution<size_t>(0, next_node - 1);
        const size_t attachment_node = attachment_distribution(random_engine);
        const size_t new_nodes = std::min(
            cycle_size_distribution(random_engine) - 1, node_count - next_node);
        if (new_nodes == 1)
        {
            add_edge(attachment_node, next_node++, 2.0);
            continue;
        }

        size_t previous_node = attachment_node;
        for (size_t i = 0; i < new_nodes; ++i)
        {
            const size_t new_node = next_node++;
            add_edge(previous_node, new_node, 1.0);
            previous_node = new_node;
        }
        add_edge(previous_node, attachment_node, 1.0);
    }

    // Hide the construction order by randomly permuting vertex labels.
    auto labels = std::vector<size_t>(node_count);
    std::iota(labels.begin(), labels.end(), 0);
    std::shuffle(labels.begin(), labels.end(), random_engine);
    auto relabeled = std::vector<std::vector<std::pair<size_t, double>>>(node_count);
    for (size_t node = 0; node < node_count; ++node)
        for (const auto &[neighbor, weight] : adjacency_lists[node])
            relabeled[labels[node]].emplace_back(labels[neighbor], weight);

    auto vertices = std::vector<size_t>(node_count + 1, 0);
    auto edges = std::vector<size_t>();
    auto weights = std::vector<double>();
    for (size_t node = 0; node < node_count; ++node)
    {
        for (const auto &[neighbor, weight] : relabeled[node])
        {
            edges.push_back(neighbor);
            weights.push_back(weight);
        }
        vertices[node + 1] = edges.size();
    }
    return WeightedCRFGraph<>{{vertices, edges}, weights};
}
