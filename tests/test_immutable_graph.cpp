#include <gtest/gtest.h>

#include <filesystem>
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
