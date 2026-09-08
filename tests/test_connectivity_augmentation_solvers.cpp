#include <cmath>
#include <concepts>
#include <filesystem>
#include <memory>
#include <unordered_set>

#include <gtest/gtest.h>

#include "algorithms/global_mincut/viecut.h"
#include "data_structure/mutable_graph.h"

#undef DEBUG

#include "HeiConnect/data_structures/graph_utils.hpp"
#include "HeiConnect/conn_aug/algorithms/wheel_con.hpp"
#include "HeiConnect/conn_aug/reducers/full_single_dom_reducer.hpp"
#include "HeiConnect/sc_reduction/transform_single_builders.hpp"
#include "HeiConnect/set_cover/solver_greedy.hpp"
#include "HeiConnect/set_cover/solver_greedy_cheapest.hpp"
#include "HeiConnect/set_cover/solver_ilp.hpp"
#include "HeiConnect/set_cover/solver_lazy_block_tree_ilp.hpp"
#include "HeiConnect/set_cover/trimmer.hpp"

namespace
{
    using SelectedLinks = std::unordered_set<size_t>;

    mutableGraphPtr make_viecut_graph(
        const WeightedCRFGraph<>& graph,
        const WeightedCRFGraph<>& link_graph,
        const SelectedLinks& selected_links)
    {
        auto viecut_graph = std::make_shared<mutable_graph>();
        viecut_graph->start_construction(graph.num_vertices());
        for (size_t u = 0; u < graph.num_vertices(); ++u)
        {
            viecut_graph->new_node();
            for (size_t e = graph.graph.vertices[u]; e < graph.graph.vertices[u + 1]; ++e)
            {
                const size_t v = graph.graph.edges[e];
                if (u < v)
                {
                    viecut_graph->new_edge(u, v, static_cast<EdgeWeight>(std::llround(graph.weights[e])));
                }
            }
        }

        const auto links = link_graph.csr_to_vec_links();
        for (const size_t link_id : selected_links)
        {
            auto [u, v, cost] = links.at(link_id);
            (void)cost;
            if (u > v)
            {
                std::swap(u, v);
            }
            viecut_graph->new_edge(u, v, 1);
        }
        viecut_graph->finish_construction();
        return viecut_graph;
    }

    EdgeWeight viecut_min_cut(
        const WeightedCRFGraph<>& graph,
        const WeightedCRFGraph<>& link_graph,
        const SelectedLinks& selected_links = {})
    {
        viecut<mutableGraphPtr> min_cut_solver;
        return min_cut_solver.perform_minimum_cut(make_viecut_graph(graph, link_graph, selected_links));
    }

    template<typename SetCoverType, typename SolverType>
    SelectedLinks solve_with_basic_context(const SetCoverType& set_cover, SolverType solver)
    {
        BoundContext context{BasicContext(set_cover), USSolution{}};
        if constexpr (std::same_as<decltype(solver.solve(set_cover, context)), bool>)
        {
            EXPECT_TRUE(solver.solve(set_cover, context));
        }
        else
        {
            solver.solve(set_cover, context);
        }
        return context.get_solution();
    }

    template<typename SolverType>
    void expect_csr_solver_increases_connectivity(SolverType solver)
    {
        const auto graph = create_cycle_graph_undirected(4);
        const auto link_graph = graph.generate_links([](size_t, size_t) { return 1.0; });
        const auto set_cover = construct_set_cover(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        const auto solution = solve_with_basic_context(set_cover, std::move(solver));

        ASSERT_FALSE(solution.empty());
        EXPECT_GT(viecut_min_cut(graph, link_graph, solution), viecut_min_cut(graph, link_graph));
    }

    void expect_reconstructed_reduced_solution_increases_connectivity(
        const std::filesystem::path& graph_file,
        const std::filesystem::path& link_file,
        ConnectivityAugmentationReductionConfig config = {})
    {
        const auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
        const auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

        auto link_remap = make_identity_link_remap(link_graph);
        UnionFind uf(graph.num_vertices());
        USSolution forced_solution;
        FullSingleDomReducer<> reducer{config};
        auto [reduced_graph, reduced_link_graph, final_link_remap, final_uf] =
            reducer.run(graph, link_graph, link_remap, uf, forced_solution);
        (void)final_uf;

        const auto reduced_set_cover = construct_set_cover(
            reduced_graph.graph.vertices,
            reduced_graph.graph.edges,
            reduced_graph.weights,
            reduced_link_graph.graph.vertices,
            reduced_link_graph.graph.edges,
            reduced_link_graph.weights);
        const auto reduced_solution = solve_with_basic_context(reduced_set_cover, GreedySetCoverSolver<>{});

        SelectedLinks original_solution = forced_solution.get_solution();
        const auto reduced_links = reduced_link_graph.csr_to_vec_links();
        for (const size_t reduced_link_id : reduced_solution)
        {
            const auto& [u, v, weight] = reduced_links[reduced_link_id];
            (void)weight;
            original_solution.insert(final_link_remap.at(normalize_link(u, v)).original_id);
        }

        ASSERT_FALSE(original_solution.empty());
        EXPECT_GT(
            viecut_min_cut(graph, link_graph, original_solution),
            viecut_min_cut(graph, link_graph));
    }
}

TEST(ConnectivityAugmentationSolvers, CSRGreedyIncreasesConnectivity)
{
    expect_csr_solver_increases_connectivity(GreedySetCoverSolver<>{});
}

TEST(ConnectivityAugmentationSolvers, WheelConIncreasesCycleConnectivity)
{
    const auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = graph.generate_links([](size_t, size_t) { return 1.0; });
    auto links = std::vector<HeiConnect::conn_aug::Link<>>{};
    for (const auto& [u, v, cost] : link_graph.csr_to_vec_links())
    {
        links.push_back({u, v, cost});
    }
    using Instance = HeiConnect::conn_aug::Instance<WeightedCRFGraph<>>;
    const Instance instance{graph, std::move(links)};
    const auto solution = HeiConnect::conn_aug::WheelCon{}.solve(instance);
    const SelectedLinks selected_links(
        solution.selected_link_ids.begin(),
        solution.selected_link_ids.end());

    ASSERT_EQ(solution.status, HeiConnect::conn_aug::SolveStatus::Feasible);
    ASSERT_EQ(selected_links.size(), 2);
    EXPECT_GT(
        viecut_min_cut(graph, link_graph, selected_links),
        viecut_min_cut(graph, link_graph));
}

TEST(ConnectivityAugmentationSolvers, WheelConIncreasesTreeConnectivity)
{
    const auto graph = create_star_graph<size_t, size_t, double>(4);
    const auto link_graph = graph.generate_links([](size_t, size_t) { return 1.0; });
    auto links = std::vector<HeiConnect::conn_aug::Link<>>{};
    for (const auto& [u, v, cost] : link_graph.csr_to_vec_links())
    {
        links.push_back({u, v, cost});
    }
    using Instance = HeiConnect::conn_aug::Instance<WeightedCRFGraph<>>;
    const Instance instance{graph, std::move(links)};
    const auto solution = HeiConnect::conn_aug::WheelCon{}.solve(instance);
    const SelectedLinks selected_links(
        solution.selected_link_ids.begin(),
        solution.selected_link_ids.end());

    ASSERT_EQ(solution.status, HeiConnect::conn_aug::SolveStatus::Feasible);
    ASSERT_EQ(selected_links.size(), 2);
    EXPECT_GT(
        viecut_min_cut(graph, link_graph, selected_links),
        viecut_min_cut(graph, link_graph));
}

TEST(ConnectivityAugmentationSolvers, CSRGreedyCheapestIncreasesConnectivity)
{
    expect_csr_solver_increases_connectivity(SetCoverSolverGreedyCheapest{});
}

TEST(ConnectivityAugmentationSolvers, CSRILPIncreasesConnectivity)
{
    try
    {
        expect_csr_solver_increases_connectivity(SetCoverSolverILP<>{});
    }
    catch (const GRBException& error)
    {
        if (error.getErrorCode() == GRB_ERROR_NO_LICENSE)
        {
            GTEST_SKIP() << "Gurobi license unavailable: " << error.getMessage();
        }
        throw;
    }
}

TEST(ConnectivityAugmentationSolvers, LazyBlockTreeILPSkipsFullyReducedGraph)
{
    const auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr({{0, 2, 2.0}, {1, 3, 3.0}}, 4);
    FullSingleDomReducer<> reducer;
    UnionFind uf(graph.num_vertices());
    ConnAugLinkRemap link_remap;
    USSolution forced_solution;
    auto [reduced_graph, reduced_links, remap, reduced_uf] =
        reducer.run(graph, link_graph, link_remap, uf, forced_solution);
    ASSERT_EQ(reduced_graph.num_vertices(), 1);
    ASSERT_EQ(forced_solution.get_solution_size(), 2);

    USSolution solution;
    HeiConnect::LazyBlockTreeSolverILP<1> solver;
    ASSERT_TRUE(solver.solve(reduced_graph, reduced_links, solution));
    EXPECT_TRUE(solver.is_feasible());
    EXPECT_TRUE(solution.get_solution().empty());
    ASSERT_TRUE(solver.emit_metrics().has_value());
    EXPECT_GT(viecut_min_cut(graph, link_graph, forced_solution.get_solution()), viecut_min_cut(graph, link_graph));
}

TEST(ConnectivityAugmentationSolvers, ZeroCutSetCoverNeedsNoSolverWork)
{
    const WeightedCRFGraph<> graph{{{0, 0}, {}}, {}};
    const auto set_cover = construct_set_cover(
        graph.graph.vertices, graph.graph.edges, graph.weights,
        graph.graph.vertices, graph.graph.edges, graph.weights);
    EXPECT_EQ(set_cover.get_num_elements(), 0);
    EXPECT_EQ(set_cover.get_num_sets(), 0);
    EXPECT_TRUE(solve_with_basic_context(set_cover, GreedySetCoverSolver<>{}).empty());
    EXPECT_TRUE(solve_with_basic_context(set_cover, SetCoverSolverGreedyCheapest{}).empty());
    EXPECT_TRUE(solve_with_basic_context(set_cover, SetCoverSolverILP<1>{}).empty());

    const auto bit_set_cover = construct_set_cover_bit_matrix(
        graph.graph.vertices, graph.graph.edges, graph.weights,
        graph.graph.vertices, graph.graph.edges, graph.weights);
    EXPECT_EQ(bit_set_cover.get_num_elements(), 0);
    EXPECT_EQ(bit_set_cover.get_num_sets(), 0);
    EXPECT_EQ(bit_set_cover.get_n_cols(), 0);
}

TEST(ConnectivityAugmentationSolvers, LazyBlockTreeILPIncreasesConnectivity)
{
    try
    {
        const auto graph = create_cycle_graph_undirected(4);
        const auto link_graph = graph.generate_links([](size_t, size_t) { return 1.0; });
        USSolution solution;
        HeiConnect::LazyBlockTreeSolverILP<> solver;

        ASSERT_TRUE(solver.solve(graph, link_graph, solution));
        ASSERT_FALSE(solution.get_solution().empty());
        EXPECT_GT(
            viecut_min_cut(graph, link_graph, solution.get_solution()),
            viecut_min_cut(graph, link_graph));
    }
    catch (const GRBException& error)
    {
        if (error.getErrorCode() == GRB_ERROR_NO_LICENSE)
        {
            GTEST_SKIP() << "Gurobi license unavailable: " << error.getMessage();
        }
        throw;
    }
}

TEST(ConnectivityAugmentationSolvers, ReconstructedReducedCycleSolutionIncreasesOriginalConnectivity)
{
    const auto project_dir = std::filesystem::path(__FILE__).parent_path().parent_path();
    expect_reconstructed_reduced_solution_increases_connectivity(
        project_dir / "datasets/cycles/cycle_50.xml",
        project_dir / "datasets/cycles/cycle_50-float_uniform_0_1.links");
}

TEST(ConnectivityAugmentationSolvers, ReconstructedReducedCycleSolutionWithoutElementDomination)
{
    const auto project_dir = std::filesystem::path(__FILE__).parent_path().parent_path();
    ConnectivityAugmentationReductionConfig config;
    config.run_element_domination = false;
    expect_reconstructed_reduced_solution_increases_connectivity(
        project_dir / "datasets/cycles/cycle_50.xml",
        project_dir / "datasets/cycles/cycle_50-float_uniform_0_1.links",
        config);
}

TEST(ConnectivityAugmentationSolvers, ReconstructedReducedCycleSolutionWithoutSingleLink)
{
    const auto project_dir = std::filesystem::path(__FILE__).parent_path().parent_path();
    ConnectivityAugmentationReductionConfig config;
    config.run_single_link = false;
    expect_reconstructed_reduced_solution_increases_connectivity(
        project_dir / "datasets/cycles/cycle_50.xml",
        project_dir / "datasets/cycles/cycle_50-float_uniform_0_1.links",
        config);
}

TEST(ConnectivityAugmentationSolvers, ReconstructedReducedStarSolutionIncreasesOriginalConnectivity)
{
    const auto project_dir = std::filesystem::path(__FILE__).parent_path().parent_path();
    expect_reconstructed_reduced_solution_increases_connectivity(
        project_dir / "datasets/stars/star_50.xml",
        project_dir / "datasets/stars/star_50-float_uniform_0_1.links");
}

TEST(ConnectivityAugmentationSolvers, ReconstructedReducedTreeSolutionIncreasesOriginalConnectivity)
{
    const auto project_dir = std::filesystem::path(__FILE__).parent_path().parent_path();
    expect_reconstructed_reduced_solution_increases_connectivity(
        project_dir / "datasets/trees/tree_50.xml",
        project_dir / "datasets/trees/tree_50.links");
}

TEST(ConnectivityAugmentationSolvers, OracleGreedyIncreasesConnectivity)
{
    const auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = graph.generate_links([](size_t, size_t) { return 1.0; });
    const auto set_cover = construct_set_cover_oracle(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    const auto solution = solve_with_basic_context(set_cover, GreedySetCoverSolver<>{});

    ASSERT_FALSE(solution.empty());
    EXPECT_GT(viecut_min_cut(graph, link_graph, solution), viecut_min_cut(graph, link_graph));
}

TEST(ConnectivityAugmentationSolvers, DoubleGreedyIncreasesConnectivity)
{
    const auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = graph.generate_links([](size_t, size_t) { return 1.0; });
    const auto set_cover = construct_set_cover_cyc_pseudo_ancestry_vec(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    BitPackedContext pseudo_context{set_cover.first()};
    CycContext cycle_context{set_cover.second()};
    ContextDouble double_context{set_cover, std::move(pseudo_context), std::move(cycle_context)};
    BoundContext context{std::move(double_context), USSolution{}};
    GreedySetCoverSolver<> solver;
    solver.solve(set_cover, context);

    ASSERT_FALSE(context.get_solution().empty());
    EXPECT_GT(
        viecut_min_cut(graph, link_graph, context.get_solution()),
        viecut_min_cut(graph, link_graph));
}

TEST(ConnectivityAugmentationSolvers, DoubleGreedyCheapestIncreasesConnectivity)
{
    const auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = graph.generate_links([](size_t, size_t) { return 1.0; });
    const auto set_cover = construct_set_cover_cyc_pseudo_ancestry_vec(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    BitPackedContext pseudo_context{set_cover.first()};
    CycContext cycle_context{set_cover.second()};
    ContextDouble double_context{set_cover, std::move(pseudo_context), std::move(cycle_context)};
    BoundContext context{std::move(double_context), USSolution{}};
    SetCoverSolverGreedyCheapest solver;
    solver.solve(set_cover, context);

    ASSERT_FALSE(context.get_solution().empty());
    EXPECT_EQ(context.get_total_covered_elements(), set_cover.get_num_elements());
    EXPECT_GT(
        viecut_min_cut(graph, link_graph, context.get_solution()),
        viecut_min_cut(graph, link_graph));
}

TEST(ConnectivityAugmentationSolvers, DoubleTrimmerReplaysRetainedSolution)
{
    const auto graph = create_cycle_graph_undirected(6);
    const auto link_graph = graph.generate_links([](size_t, size_t) { return 1.0; });
    auto set_cover = construct_set_cover_cyc_pseudo_ancestry_vec(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    auto set_cover_ptr = std::make_shared<const decltype(set_cover)>(std::move(set_cover));
    BitPackedContext pseudo_context{set_cover_ptr->first()};
    CycContext cycle_context{set_cover_ptr->second()};
    ContextDouble double_context{*set_cover_ptr, std::move(pseudo_context), std::move(cycle_context)};
    BoundContext context{std::move(double_context), USSolution{}};
    for (size_t set_index = 0; set_index < set_cover_ptr->get_num_sets(); ++set_index)
    {
        context.add_set(set_index);
    }
    const size_t untrimmed_size = context.get_solution_size();
    ASSERT_EQ(context.get_total_covered_elements(), set_cover_ptr->get_num_elements());

    SetCoverTrimmer trimmer;
    auto [trimmed_set_cover, trimmed_context] = trimmer(set_cover_ptr, std::move(context));

    EXPECT_LT(trimmed_context.get_solution_size(), untrimmed_size);
    EXPECT_EQ(trimmed_context.get_total_covered_elements(), trimmed_set_cover->get_num_elements());
    EXPECT_GT(
        viecut_min_cut(graph, link_graph, trimmed_context.get_solution()),
        viecut_min_cut(graph, link_graph));
}
