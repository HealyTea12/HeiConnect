#include <gtest/gtest.h>

#include <filesystem>
#include <limits>
#include <set>

#include "HeiConnect/data_structures/graph_utils.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"

TEST(AddLinks, General)
{
    CRFGraph graph{{0, 2, 4, 4}, {1, 2, 0, 2}}; // edges: 0->1,0->2,1->0,1->2
    std::vector<double> weights = {1.0, 1.0, 1.0, 1.0};
    WeightedCRFGraph wgraph{graph, weights};

    CRFGraph link_graph{{0, 0, 0, 2}, {0, 1}}; // edges: 2->0, 2->1
    std::vector<double> link_weights = {2.0, 2.0, 2.0};
    WeightedCRFGraph wlink_graph{link_graph, link_weights};

    std::unordered_set<size_t> selected_edges = {0, 1};

    auto new_graph = wgraph.add_links(wlink_graph, selected_edges);

    std::vector<size_t> expected_vertices = {0, 2, 4, 6};
    std::vector<size_t> expected_edges = {1, 2, 0, 2, 0, 1};
    std::vector<double> expected_weights = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0};

    EXPECT_EQ(new_graph.graph.vertices, expected_vertices);
    EXPECT_EQ(new_graph.graph.edges, expected_edges);
    EXPECT_EQ(new_graph.weights, expected_weights);
}

TEST(ReadFromFileGraphML, General)
{
    auto src_dir = std::filesystem::current_path().parent_path();
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(src_dir / "tests/data/k4.xml");
    ASSERT_EQ(graph.graph.vertices.size(), 5);
    ASSERT_EQ(graph.graph.edges.size(), 12);
    ASSERT_EQ(graph.weights.size(), 12);
    auto expected_vertices = std::vector<size_t>{0, 3, 6, 9, 12};
    EXPECT_EQ(graph.graph.vertices, expected_vertices);
    auto expected_edges = std::vector<size_t>{
        1, 2, 3,
        0, 2, 3,
        0, 1, 3,
        0, 1, 2};
    EXPECT_EQ(graph.graph.edges, expected_edges);
    auto expected_weights = std::vector<double>{
        1.0, 1.0, 1.0,
        1.0, 1.0, 1.0,
        1.0, 1.0, 1.0,
        1.0, 1.0, 1.0};
    EXPECT_EQ(graph.weights, expected_weights);
}

TEST(ImmutableGraphBlockTree, Cycle5BecomesStar)
{
    auto cycle = create_cycle_graph_undirected(5);
    auto [block_tree, cycle_positions] = cycle.cactus_generate_block_tree(0);

    ASSERT_EQ(block_tree.num_vertices(), 6);
    ASSERT_EQ(cycle_positions.size(), 1);
    ASSERT_EQ(cycle_positions[0].size(), 5);

    const size_t cycle_node = 5;

    for (size_t u = 0; u < 5; ++u)
    {
        std::set<size_t> neighbors;
        for (size_t e = block_tree.graph.vertices[u]; e < block_tree.graph.vertices[u + 1]; ++e)
        {
            neighbors.insert(block_tree.graph.edges[e]);
        }
        EXPECT_EQ(neighbors.size(), 1);
        EXPECT_EQ(*neighbors.begin(), cycle_node);
    }

    std::set<size_t> center_neighbors;
    for (size_t e = block_tree.graph.vertices[cycle_node]; e < block_tree.graph.vertices[cycle_node + 1]; ++e)
    {
        center_neighbors.insert(block_tree.graph.edges[e]);
    }
    EXPECT_EQ(center_neighbors.size(), 5);
    for (size_t u = 0; u < 5; ++u)
    {
        EXPECT_TRUE(center_neighbors.contains(u));
    }
}

TEST(ImmutableGraphCactusTraversal, TraversesBridgesTwiceAndCycleEdgesOnce)
{
    const auto cactus = WeightedCRFGraph<>{
        {std::vector<size_t>{0, 2, 4, 7, 8},
         std::vector<size_t>{1, 2, 0, 2, 0, 1, 3, 2}},
        std::vector<double>(8, 1.0)};
    auto order = std::vector<size_t>{};

    cactus.cactus_for_each_hamiltonian_vertex(
        0,
        [&](size_t u) { order.push_back(u); });

    ASSERT_EQ(order.size(), 5);
    EXPECT_EQ(std::count(order.begin(), order.end(), 0), 1);
    EXPECT_EQ(std::count(order.begin(), order.end(), 1), 1);
    EXPECT_EQ(std::count(order.begin(), order.end(), 2), 2);
    EXPECT_EQ(std::count(order.begin(), order.end(), 3), 1);
    for (size_t i = 0; i < order.size(); ++i)
    {
        EXPECT_TRUE(cactus.is_edge(order[i], order[(i + 1) % order.size()]));
    }
}

// Verifies the random cactus generator's CSR, edge-count, weight, and connectivity invariants.
TEST(RandomCactus, GeneratesValidWeightedCsrGraph)
{
    constexpr size_t node_count = 25;
    constexpr size_t cycle_count = 4;
    constexpr size_t cycle_length = 4;
    const auto cactus = create_random_cactus(node_count, cycle_count, cycle_length, 42);

    EXPECT_EQ(cactus.num_vertices(), node_count);
    EXPECT_EQ(cactus.num_edges() / 2, node_count - 1 + cycle_count);
    EXPECT_EQ(cactus.graph.vertices.front(), 0);
    EXPECT_EQ(cactus.graph.vertices.back(), cactus.graph.edges.size());
    EXPECT_TRUE(std::is_sorted(cactus.graph.vertices.begin(), cactus.graph.vertices.end()));

    size_t directed_cycle_edge_count = 0;
    size_t directed_bridge_edge_count = 0;
    for (size_t edge = 0; edge < cactus.graph.edges.size(); ++edge)
    {
        EXPECT_LT(cactus.graph.edges[edge], node_count);
        if (cactus.weights[edge] == 1.0)
            ++directed_cycle_edge_count;
        if (cactus.weights[edge] == 2.0)
            ++directed_bridge_edge_count;
    }
    EXPECT_EQ(directed_cycle_edge_count, 2 * cycle_count * cycle_length);
    EXPECT_EQ(
        directed_bridge_edge_count,
        2 * (node_count - 1 - cycle_count * (cycle_length - 1)));

    const auto [block_tree, cycle_positions] = cactus.cactus_generate_block_tree(0);
    EXPECT_EQ(block_tree.num_vertices(), node_count + cycle_count);
    ASSERT_EQ(cycle_positions.size(), cycle_count);
    for (const auto &positions : cycle_positions)
    {
        EXPECT_EQ(
            std::count_if(positions.begin(), positions.end(), [](int position)
                          { return position >= 0; }),
            cycle_length);
    }

    std::vector<bool> visited(node_count, false);
    std::vector<size_t> nodes_to_visit = {0};
    visited[0] = true;
    while (!nodes_to_visit.empty())
    {
        const size_t node = nodes_to_visit.back();
        nodes_to_visit.pop_back();
        for (size_t edge = cactus.graph.vertices[node];
             edge < cactus.graph.vertices[node + 1];
             ++edge)
        {
            const size_t neighbor = cactus.graph.edges[edge];
            if (!visited[neighbor])
            {
                visited[neighbor] = true;
                nodes_to_visit.push_back(neighbor);
            }
        }
    }
    EXPECT_TRUE(std::all_of(visited.begin(), visited.end(), [](bool was_visited)
                            { return was_visited; }));
}

// Verifies that invalid cactus parameters are rejected by the CSR generator itself.
TEST(RandomCactus, RejectsInvalidParameters)
{
    EXPECT_THROW(create_random_cactus(0, 0, 3), std::invalid_argument);
    EXPECT_THROW(create_random_cactus(5, 1, 2), std::invalid_argument);
    EXPECT_THROW(create_random_cactus(5, 3, 3), std::invalid_argument);
}

TEST(ExactCycleCountCactus, RespectsCountsStructureAndWeights)
{
    for (size_t node_count = 1; node_count <= 40; ++node_count)
    {
        for (size_t cycle_count = node_count > 1 ? 1 : 0; cycle_count < node_count; ++cycle_count)
        {
            for (unsigned int seed = 0; seed < 3; ++seed)
            {
                SCOPED_TRACE(::testing::Message() << "nodes=" << node_count
                             << " cycles=" << cycle_count << " seed=" << seed);
                const auto cactus = create_random_cactus_with_cycle_count(node_count, cycle_count, seed);
                ASSERT_EQ(cactus.num_vertices(), node_count);
                const size_t bridges = std::count(cactus.weights.begin(), cactus.weights.end(), 2.0) / 2;
                ASSERT_LE(bridges, cycle_count);
                ASSERT_EQ(cactus.num_edges(), 2 * (node_count - 1 + cycle_count - bridges));
                ASSERT_EQ(cactus.weights.size(), cactus.num_edges());
                EXPECT_EQ(cactus.graph.vertices.front(), 0);
                EXPECT_EQ(cactus.graph.vertices.back(), cactus.num_edges());
                EXPECT_TRUE(std::is_sorted(cactus.graph.vertices.begin(), cactus.graph.vertices.end()));

                const auto [parent, depth] = cactus.graph.rooted_parent_depth();
                EXPECT_EQ(std::count(parent.begin(), parent.end(), node_count), 0);
                const auto [block_tree, cycles] = cactus.cactus_generate_block_tree(0);
                ASSERT_EQ(cycles.size() + bridges, cycle_count);
                EXPECT_EQ(block_tree.num_vertices(), node_count + cycles.size());
                EXPECT_EQ(block_tree.num_edges() / 2, block_tree.num_vertices() - 1);

                std::set<std::pair<size_t, size_t>> cycle_edges;
                for (const auto &positions : cycles)
                {
                    std::vector<size_t> ordered_cycle(node_count);
                    size_t cycle_size = 0;
                    for (size_t node = 0; node < node_count; ++node)
                    {
                        if (positions[node] >= 0)
                        {
                            ASSERT_LT(static_cast<size_t>(positions[node]), node_count);
                            ordered_cycle[positions[node]] = node;
                            ++cycle_size;
                        }
                    }
                    ASSERT_GE(cycle_size, 3);
                    for (size_t i = 0; i < cycle_size; ++i)
                    {
                        const size_t u = ordered_cycle[i];
                        const size_t v = ordered_cycle[(i + 1) % cycle_size];
                        EXPECT_TRUE(cactus.is_edge(u, v));
                        EXPECT_TRUE(cycle_edges.emplace(u, v).second);
                        EXPECT_TRUE(cycle_edges.emplace(v, u).second);
                    }
                }
                for (size_t u = 0; u < node_count; ++u)
                {
                    std::set<size_t> neighbors;
                    for (size_t edge = cactus.graph.vertices[u]; edge < cactus.graph.vertices[u + 1]; ++edge)
                    {
                        const size_t v = cactus.graph.edges[edge];
                        ASSERT_LT(v, node_count);
                        EXPECT_NE(u, v);
                        EXPECT_TRUE(neighbors.insert(v).second);
                        EXPECT_TRUE(cactus.is_edge(v, u));
                        EXPECT_EQ(cactus.weights[edge], cycle_edges.contains({u, v}) ? 1.0 : 2.0);
                    }
                }
            }
        }
    }
}

TEST(ExactCycleCountCactus, ReproducesSeedsAndVariesCycleSizes)
{
    const auto first = create_random_cactus_with_cycle_count(100, 60, 42);
    const auto repeated = create_random_cactus_with_cycle_count(100, 60, 42);
    const auto other = create_random_cactus_with_cycle_count(100, 60, 43);
    EXPECT_EQ(first.graph.vertices, repeated.graph.vertices);
    EXPECT_EQ(first.graph.edges, repeated.graph.edges);
    EXPECT_EQ(first.weights, repeated.weights);
    EXPECT_NE(first.graph.edges, other.graph.edges);
    const auto [block_tree, cycles] = first.cactus_generate_block_tree(0);
    std::set<size_t> sizes;
    for (const auto &positions : cycles)
        sizes.insert(std::count_if(positions.begin(), positions.end(), [](int position)
                                  { return position >= 0; }));
    EXPECT_GT(sizes.size(), 1);
    EXPECT_GT(std::count(first.weights.begin(), first.weights.end(), 2.0), 0);
}

TEST(ExactCycleCountCactus, RejectsInvalidParameters)
{
    EXPECT_THROW(create_random_cactus_with_cycle_count(0, 0), std::invalid_argument);
    EXPECT_THROW(create_random_cactus_with_cycle_count(1, 1), std::invalid_argument);
    EXPECT_THROW(create_random_cactus_with_cycle_count(2, 0), std::invalid_argument);
    EXPECT_THROW(create_random_cactus_with_cycle_count(10, 10), std::invalid_argument);
    EXPECT_THROW(create_random_cactus_with_cycle_count(10, std::numeric_limits<size_t>::max()),
                 std::invalid_argument);
}

TEST(ExactCycleCountCactus, CountsEachTreeEdgeAsOneCycle)
{
    const auto tree = create_random_cactus_with_cycle_count(10, 9);
    EXPECT_EQ(tree.num_vertices(), 10);
    EXPECT_EQ(tree.num_edges(), 18);
    EXPECT_EQ(tree.weights, std::vector<double>(18, 2.0));
    const auto [tree_blocks, tree_cycles] = tree.cactus_generate_block_tree(0);
    EXPECT_TRUE(tree_cycles.empty());

    const auto cycle = create_random_cactus_with_cycle_count(10, 1);
    EXPECT_EQ(cycle.num_vertices(), 10);
    EXPECT_EQ(cycle.num_edges(), 20);
    EXPECT_EQ(cycle.weights, std::vector<double>(20, 1.0));
    const auto [cycle_blocks, cycles] = cycle.cactus_generate_block_tree(0);
    EXPECT_EQ(cycles.size(), 1);
}

TEST(VariableCycleCactus, RespectsBudgetAndCactusStructure)
{
    for (size_t node_count = 1; node_count <= 60; ++node_count)
    {
        for (unsigned int seed = 0; seed < 5; ++seed)
        {
            SCOPED_TRACE(::testing::Message() << "nodes=" << node_count << " seed=" << seed);
            const auto cactus = create_random_cactus_with_cycle_sizes(node_count, 5, 9, seed);
            ASSERT_EQ(cactus.num_vertices(), node_count);
            ASSERT_EQ(cactus.weights.size(), cactus.num_edges());
            EXPECT_EQ(cactus.graph.vertices.front(), 0);
            EXPECT_EQ(cactus.graph.vertices.back(), cactus.num_edges());
            EXPECT_TRUE(std::is_sorted(cactus.graph.vertices.begin(), cactus.graph.vertices.end()));
            for (size_t u = 0; u < node_count; ++u)
            {
                std::set<size_t> neighbors;
                for (size_t edge = cactus.graph.vertices[u]; edge < cactus.graph.vertices[u + 1]; ++edge)
                {
                    const size_t v = cactus.graph.edges[edge];
                    ASSERT_LT(v, node_count);
                    EXPECT_NE(u, v);
                    EXPECT_TRUE(neighbors.insert(v).second);
                    EXPECT_TRUE(cactus.is_edge(v, u));
                    EXPECT_TRUE(cactus.weights[edge] == 1.0 || cactus.weights[edge] == 2.0);
                }
            }
            const auto [parent, depth] = cactus.graph.rooted_parent_depth();
            EXPECT_EQ(std::count(parent.begin(), parent.end(), node_count), 0);

            const auto [block_tree, positions] = cactus.cactus_generate_block_tree(0);
            EXPECT_EQ(block_tree.num_edges() / 2, block_tree.num_vertices() - 1);
            EXPECT_EQ(block_tree.num_vertices(), node_count + positions.size());
            size_t cycle_edges = 0;
            size_t short_cycles = 0;
            for (const auto &cycle : positions)
            {
                const size_t size = std::count_if(cycle.begin(), cycle.end(), [](int position)
                                                 { return position >= 0; });
                EXPECT_GE(size, 3);
                EXPECT_LE(size, 9);
                short_cycles += size < 5;
                cycle_edges += size;
            }
            EXPECT_LE(short_cycles, 1);
            const size_t directed_bridges = std::count(cactus.weights.begin(), cactus.weights.end(), 2.0);
            EXPECT_TRUE(directed_bridges == 0 || directed_bridges == 2);
            EXPECT_EQ(cactus.num_edges(), 2 * cycle_edges + directed_bridges);
            EXPECT_EQ(cactus.num_edges() / 2, node_count - 1 + positions.size());
        }
    }
}

TEST(VariableCycleCactus, TruncatesLastCycleAndHandlesSingleRemainingVertex)
{
    const auto truncated = create_random_cactus_with_cycle_sizes(8, 5, 5);
    const auto [block_tree, positions] = truncated.cactus_generate_block_tree(0);
    std::multiset<size_t> sizes;
    for (const auto &cycle : positions)
        sizes.insert(std::count_if(cycle.begin(), cycle.end(), [](int position)
                                   { return position >= 0; }));
    EXPECT_EQ(sizes, (std::multiset<size_t>{4, 5}));
    EXPECT_EQ(std::count(truncated.weights.begin(), truncated.weights.end(), 2.0), 0);

    const auto with_bridge = create_random_cactus_with_cycle_sizes(6, 5, 5);
    EXPECT_EQ(with_bridge.num_vertices(), 6);
    EXPECT_EQ(std::count(with_bridge.weights.begin(), with_bridge.weights.end(), 2.0), 2);

    const auto tiny = create_random_cactus_with_cycle_sizes(2, 5, 9);
    EXPECT_EQ(tiny.num_vertices(), 2);
    EXPECT_EQ(tiny.weights, (std::vector<double>{2.0, 2.0}));
}

TEST(VariableCycleCactus, ReproducesSeedsAndProducesDifferentCycleSizes)
{
    const auto first = create_random_cactus_with_cycle_sizes(400, 3, 7, 42);
    const auto repeated = create_random_cactus_with_cycle_sizes(400, 3, 7, 42);
    const auto other = create_random_cactus_with_cycle_sizes(400, 3, 7, 43);
    EXPECT_EQ(first.graph.vertices, repeated.graph.vertices);
    EXPECT_EQ(first.graph.edges, repeated.graph.edges);
    EXPECT_EQ(first.weights, repeated.weights);
    EXPECT_NE(first.graph.edges, other.graph.edges);
    const auto [block_tree, positions] = first.cactus_generate_block_tree(0);
    std::set<size_t> sizes;
    for (const auto &cycle : positions)
        sizes.insert(std::count_if(cycle.begin(), cycle.end(), [](int position)
                                   { return position >= 0; }));
    EXPECT_EQ(sizes, (std::set<size_t>{3, 4, 5, 6, 7}));
}

TEST(VariableCycleCactus, RejectsInvalidParameters)
{
    EXPECT_THROW(create_random_cactus_with_cycle_sizes(0, 3, 5), std::invalid_argument);
    EXPECT_THROW(create_random_cactus_with_cycle_sizes(10, 1, 5), std::invalid_argument);
    EXPECT_THROW(create_random_cactus_with_cycle_sizes(10, 5, 3), std::invalid_argument);
}

TEST(VariableCycleCactus, SizeTwoProducesBridges)
{
    const auto tree = create_random_cactus_with_cycle_sizes(40, 2, 2);
    EXPECT_EQ(tree.num_vertices(), 40);
    EXPECT_EQ(tree.num_edges(), 2 * 39);
    EXPECT_TRUE(std::all_of(tree.weights.begin(), tree.weights.end(), [](double weight)
                            { return weight == 2.0; }));
    const auto [parent, depth] = tree.graph.rooted_parent_depth();
    EXPECT_EQ(std::count(parent.begin(), parent.end(), 40), 0);
    const auto [tree_blocks, tree_cycles] = tree.cactus_generate_block_tree(0);
    EXPECT_TRUE(tree_cycles.empty());

    const auto mixed = create_random_cactus_with_cycle_sizes(400, 2, 7);
    EXPECT_EQ(mixed.num_vertices(), 400);
    EXPECT_GT(std::count(mixed.weights.begin(), mixed.weights.end(), 2.0), 2);
    EXPECT_GT(std::count(mixed.weights.begin(), mixed.weights.end(), 1.0), 0);
    const auto [block_tree, cycles] = mixed.cactus_generate_block_tree(0);
    EXPECT_EQ(block_tree.num_edges() / 2, block_tree.num_vertices() - 1);
    EXPECT_EQ(mixed.num_edges() / 2, 399 + cycles.size());
}
