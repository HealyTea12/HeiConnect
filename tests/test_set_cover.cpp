
#include <gtest/gtest.h>

#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/sc_reduction/transform_single_builders.hpp"
#include "HeiConnect/min_cut/simple_mincut.hpp"
#include "HeiConnect/viecut_runner.hpp" // Really would like not to need this
#include "HeiConnect/data_structures/graph_utils.hpp"

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
    SetCover sc{a, b, costs, 4};

    GreedySetCoverSolver<SetCover, 0> solver;
    USSolution solution;
    solver.solve(sc, solution);

    std::unordered_set<size_t> expected_solution = {1, 2};
    EXPECT_EQ(solution.get_solution(), expected_solution);
}

TEST(Transpose, Basic)
{
    std::vector<uint8_t> input = {
        0b00001111, // row 0
        0b00110011, // row 1
        0b01010101  // row 2
    };
    auto output = transpose(input, 3, 1);
    std::vector<uint8_t> expected_output = {
        0b00000111, // col 0
        0b00000011, // col 1
        0b00000101, // col 2
        0b00000001, // col 3
        0b00000110, // col 4
        0b00000010, // col 5
        0b00000100, // col 6
        0b00000000  // col 7
    };
    EXPECT_EQ(expected_output, output);
}

TEST(GraphAugmentationToSetCover, Cycle)
{
    auto cycle_graph = create_cycle_graph_undirected(10); // 0 <-> 1 <-> 2 <-> 3 ... <-> 0
    auto link_graph = cycle_graph.generate_links([](size_t u, size_t v)
                                                 { return 1.0; });

    auto sc = construct_set_cover(
        cycle_graph.graph.vertices,
        cycle_graph.graph.edges,
        cycle_graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    for (size_t u = 0; u < link_graph.graph.vertices.size() - 1; u++)
    {
        for (size_t e = link_graph.graph.vertices[u]; e < link_graph.graph.vertices[u + 1]; e++)
        {
            size_t v = link_graph.graph.edges[e];
            size_t n_cuts_covered = 0;
            sc.forEachElement(e, [&](size_t)
                              { ++n_cuts_covered; });
            auto distance_u_v = (std::max(u, v) - std::min(u, v));
            auto expected_cuts_covered = distance_u_v * ((cycle_graph.graph.vertices.size() - 1) - distance_u_v);
            ASSERT_EQ(n_cuts_covered, expected_cuts_covered);
        }
    }
}

// TEST(GraphAugmentationToSetCover, FromFile)
//{
//     auto src_dir = std::filesystem::current_path().parent_path();
//     std::filesystem::recursive_directory_iterator it(src_dir / "graphs");
//     for (const auto &entry : it)
//     {
//         if (entry.is_directory() || entry.path().extension() != ".graph")
//             continue;
//         std::cout << "Testing graph: " << entry.path() << std::endl;
//         std::stringstream command{};
//         command << src_dir / "extern" / "VieCut" / "build" / "mincut"
//                 << " " << entry.path()
//                 << " cactus -t " << entry.path().string() + ".cactus";
//         system(command.str().c_str());
//
//         auto cactus_path = entry.path().parent_path() / (entry.path().filename().string() + ".cactus");
//         auto graph = WeightedCRFGraph<>::read_from_file_graphML(cactus_path);
//         ASSERT_FALSE(graph.graph.vertices.empty());
//         ASSERT_FALSE(graph.graph.edges.empty());
//
//         auto link_graph = graph.generate_links([](size_t u, size_t v)
//                                                { return 1.0; });
//         if (link_graph.graph.edges.empty())
//         {
//             std::cout << "Graph " << entry.path() << " is complete, skipping." << std::endl;
//             continue;
//         }
//         std::cout << "Constructing set cover...\n";
//         auto set_cover_instance = construct_set_cover(
//             graph.graph.vertices,
//             graph.graph.edges,
//             graph.weights,
//             link_graph.graph.vertices,
//             link_graph.graph.edges,
//             link_graph.weights);
//
//         ASSERT_TRUE(set_cover_instance.a.size() == link_graph.graph.edges.size() + 1);
//         std::cout << "Solving set cover...\n";
//         auto solver = SetCoverSolverGreedyParallel(set_cover_instance);
//         solver.solve();
//         auto solution = solver.get_solution();
//         ASSERT_FALSE(solution.empty());
//         // calculate min_cut of original graph
//         auto original_min_cut = global_mincut_simple(graph);
//         // construct new graph with add links from solution
//         auto new_graph = graph.add_links(
//             link_graph,
//             solution);
//         new_graph = new_graph.make_bidirectional();
//         auto new_min_cut = global_mincut_simple(new_graph);
//         // check if min_cut has increased
//         ASSERT_GT(new_min_cut, original_min_cut);
//     }
// }
