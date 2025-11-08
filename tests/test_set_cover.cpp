
#include <gtest/gtest.h>
#include "set_cover/set_cover.hpp"
#include "set_cover/graph.hpp"
#include "set_cover/transform_single.hpp"
#include "set_cover/read_file.hpp"
#include "min_cut/simple_mincut.hpp"

TEST(Utils, Max)
{
    std::vector<int> a = {2, 10, 1, 55, -2, -10, 28};
    ASSERT_EQ(max(a), 55);
}

TEST(Utils, ArgMin)
{
    std::vector<double> a = {2.5, 10.1, 1.0, 55.2, -2.3, -10.4, 28.9};
    ASSERT_EQ(argmin(a), 5);
    ASSERT_EQ(argmax(a), 3);
}

TEST(SetCoverTest, SetCoverBasic)
{
    std::vector<size_t> a = {0, 2, 4, 7};
    std::vector<size_t> b = {0, 1, 1, 2, 0, 2, 3};
    std::vector<double> costs = {3.0, 2.0, 4.0};
    SetCover sc{a, b, costs};

    SetCoverSolverGreedyParallel solver(sc);
    solver.solve();
    auto solution = solver.get_solution();

    std::unordered_set<size_t> expected_solution = {1, 2};
    EXPECT_EQ(solution, expected_solution);
}

TEST(AddLinks, General)
{
    CRFGraph graph{{0, 2, 4, 4}, {1, 2, 0, 2}}; // edges: 0->1,0->2,1->0,1->2
    std::vector<double> weights = {1.0, 1.0, 1.0, 1.0};
    WeightedCRFGraph wgraph{graph, weights};

    CRFGraph link_graph{{0, 0, 0, 2}, {0, 1}}; // edges: 2->0, 2->1
    std::vector<double> link_weights = {2.0, 2.0, 2.0};
    WeightedCRFGraph wlink_graph{link_graph, link_weights};

    std::unordered_set<size_t> selected_edges = {0, 1}; // select edge 1->0

    auto new_graph = add_links(wgraph, wlink_graph, selected_edges);

    // new graph should have edges: 0->1,0->2,1->0,1->2,2->0,2->1
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
    auto graph = read_from_file_graphML(src_dir / "tests/data/k4.xml");
    ASSERT_EQ(graph.graph.vertices.size(), 5);
    ASSERT_EQ(graph.graph.edges.size(), 12);
    ASSERT_EQ(graph.weights.size(), 12);
    auto expected_vertices = std::vector<size_t>{0, 3, 6, 9, 12};
    EXPECT_EQ(graph.graph.vertices, expected_vertices);
    auto expected_edges = std::vector<size_t>{
        1, 2, 3, // edges from vertex 0
        0, 2, 3, // edges from vertex 1
        0, 1, 3, // edges from vertex 2
        0, 1, 2  // edges from vertex 3
    }; // might not be in this order though?
    EXPECT_EQ(graph.graph.edges, expected_edges);
    auto expected_weights = std::vector<double>{
        1.0, 1.0, 1.0, // weights for edges from vertex 0
        1.0, 1.0, 1.0, // weights for edges from vertex 1
        1.0, 1.0, 1.0, // weights for edges from vertex 2
        1.0, 1.0, 1.0  // weights for edges from vertex 3
    }; // might not be in this order though?
    EXPECT_EQ(graph.weights, expected_weights);
}

TEST(GraphAugmentationToSetCover, FromFile)
{
    auto src_dir = std::filesystem::current_path().parent_path();
    auto graph = read_from_file_graphML(src_dir / "graphs/dimacs10/clustering/email.cactus.xml");
    ASSERT_FALSE(graph.graph.vertices.empty());
    ASSERT_FALSE(graph.graph.edges.empty());

    auto link_graph = generate_links(graph, [](size_t u, size_t v)
                                     { return 1.0; });
    auto set_cover_instance = construct_set_cover(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);

    ASSERT_TRUE(set_cover_instance.a.size() == link_graph.graph.edges.size() + 1);
    auto solver = SetCoverSolverGreedyParallel(set_cover_instance);
    solver.solve();
    auto solution = solver.get_solution();
    ASSERT_FALSE(solution.empty());
    // calculate min_cut of original graph
    auto original_min_cut = global_mincut_simple(graph);
    // construct new graph with add links from solution
    auto new_graph = add_links(graph,
                               link_graph,
                               solution);
    new_graph = make_bidirectional(new_graph);
    auto new_min_cut = global_mincut_simple(new_graph);
    // check if min_cut has increased
    ASSERT_GT(new_min_cut, original_min_cut);
}