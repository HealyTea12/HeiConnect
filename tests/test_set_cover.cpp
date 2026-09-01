
#include <gtest/gtest.h>

#include <memory>

#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/set_cover/solver_greedy_cheapest.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/sc_reduction/transform_single_builders.hpp"
#include "HeiConnect/min_cut/simple_mincut.hpp"
#include "HeiConnect/viecut_runner.hpp" // Really would like not to need this
#include "HeiConnect/data_structures/graph_utils.hpp"

namespace
{
    using PseudoSetCover = SetCoverPseudo<size_t, size_t>;
    using CycleSetCover = SetCoverCyc<size_t, size_t, double>;
    using DoubleContextType = ContextDouble<
        PseudoSetCover,
        CycleSetCover,
        BitPackedContext<PseudoSetCover>,
        CycContext<size_t, size_t, double>>;

    static_assert(IncrementallyRemovableContext<BasicContext<SetCover<>>>);
    static_assert(!IncrementallyRemovableContext<BitPackedContext<PseudoSetCover>>);
    static_assert(!IncrementallyRemovableContext<CycContext<size_t, size_t, double>>);
    static_assert(!IncrementallyRemovableContext<DoubleContextType>);
    static_assert(IncrementallyRemovableContext<BoundContext<BasicContext<SetCover<>>, USSolution>>);
    static_assert(!IncrementallyRemovableContext<BoundContext<DoubleContextType, USSolution>>);

    struct CountingRemovalContext
    {
        bool result;
        std::shared_ptr<size_t> calls;

        bool can_remove(size_t, const USSolution&) const
        {
            ++*calls;
            return result;
        }
    };
}

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
    SetCover<> sc{a, b, costs, 4};

    GreedySetCoverSolver<0> solver;
    BoundContext context{BasicContext(sc), USSolution{}};
    solver.solve(sc, context);

    std::unordered_set<size_t> expected_solution = {1, 2};
    EXPECT_EQ(context.get_solution(), expected_solution);
}

TEST(SetCoverTest, GreedyCheapestSkipsSetsWithoutNewCoverage)
{
    std::vector<size_t> offsets = {0, 0, 1, 2};
    std::vector<size_t> elements = {0, 1};
    std::vector<double> costs = {1.0, 2.0, 3.0};
    SetCover<> set_cover{offsets, elements, costs, 2};
    BoundContext context{BasicContext(set_cover), USSolution{}};
    SetCoverSolverGreedyCheapest solver;

    solver.solve(set_cover, context);

    const std::unordered_set<size_t> expected_solution = {1, 2};
    EXPECT_EQ(context.get_solution(), expected_solution);
}

TEST(SetCoverTest, GreedyCheapestReportsInfeasibleInstance)
{
    std::vector<size_t> offsets = {0, 1};
    std::vector<size_t> elements = {0};
    std::vector<double> costs = {1.0};
    SetCover<> set_cover{offsets, elements, costs, 2};
    BoundContext context{BasicContext(set_cover), USSolution{}};
    SetCoverSolverGreedyCheapest solver;

    EXPECT_THROW(solver.solve(set_cover, context), std::runtime_error);
}

TEST(SetCoverTest, DoubleContextChecksRemovalFromLeftToRight)
{
    SetCover<> first{{0, 1}, {0}, {1.0}, 1};
    SetCover<> second{{0, 1}, {0}, {1.0}, 1};
    SetCoverDouble set_cover{std::move(first), std::move(second)};
    const auto first_calls = std::make_shared<size_t>(0);
    const auto second_calls = std::make_shared<size_t>(0);
    CountingRemovalContext first_context{false, first_calls};
    CountingRemovalContext second_context{true, second_calls};
    ContextDouble context{set_cover, std::move(first_context), std::move(second_context)};
    const USSolution solution;

    EXPECT_FALSE(context.can_remove(0, solution));
    EXPECT_EQ(*first_calls, 1);
    EXPECT_EQ(*second_calls, 0);
}

TEST(SetCoverTest, OptionalTrimmerPreservesSolutionWhenDisabled)
{
    SetCover<> set_cover{{0, 1, 2}, {0, 0}, {2.0, 1.0}, 1};
    BoundContext context{BasicContext(set_cover), USSolution{}};
    context.add_set(0);
    context.add_set(1);
    OptionalSetCoverTrimmer<1> trimmer{false};

    trimmer.trim(set_cover, context);

    EXPECT_EQ(context.get_solution_size(), 2);
    EXPECT_EQ(context.get_total_covered_elements(), 1);
    const auto metrics = trimmer.emit_metrics();
    ASSERT_TRUE(metrics.has_value());
    ASSERT_EQ(metrics->size(), 1);
    EXPECT_EQ(metrics->front().name, "skipped");
    EXPECT_EQ(metrics->front().printable_value, "true");
}

TEST(SetCoverTest, OptionalTrimmerRemovesRedundantSetsWhenEnabled)
{
    SetCover<> set_cover{{0, 1, 2}, {0, 0}, {2.0, 1.0}, 1};
    BoundContext context{BasicContext(set_cover), USSolution{}};
    context.add_set(0);
    context.add_set(1);
    OptionalSetCoverTrimmer<1> trimmer{true};

    trimmer.trim(set_cover, context);

    EXPECT_EQ(context.get_solution_size(), 1);
    EXPECT_EQ(context.get_total_covered_elements(), 1);
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
