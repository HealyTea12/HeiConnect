#include "extern/VieCut/lib/algorithms/global_mincut/algorithms.h"

#include <gtest/gtest.h>
#include "set_cover/set_cover.hpp"
#include "set_cover/graph.hpp"
#include "set_cover/transform_single.hpp"
#include "set_cover/read_file.hpp"

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

TEST(SetCover, SetCoverBasic)
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

TEST(GraphAugmentationToSetCover, FromFile)
{
    auto graph = read_from_file("graphs/dimacs10/clustering/email.graph");
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
    SetCoverSolverGreedyParallel solver{set_cover_instance};
    solver.solve();
    auto solution = solver.get_solution();
    ASSERT_FALSE(solution.empty());
    // calculate min_cut of original graph
    auto mutable_original = graph.create_mutable_graph();
    auto original_min_cut = viecut<mutable_graph>{}.perform_minimum_cut(mutable_original);
    // construct new graph with add links from solution
    add_links(mutable_original,
              link_graph,
              solution);
    auto new_min_cut = viecut<mutable_graph>{}.perform_minimum_cut(mutable_original);
    // check if min_cut has increased
    ASSERT_GT(new_min_cut, original_min_cut);
}