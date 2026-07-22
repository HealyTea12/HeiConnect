#include <cmath>
#include <concepts>
#include <memory>
#include <unordered_set>

#include <gtest/gtest.h>

#include "algorithms/global_mincut/viecut.h"
#include "data_structure/mutable_graph.h"

#undef DEBUG

#include "HeiConnect/data_structures/graph_utils.hpp"
#include "HeiConnect/sc_reduction/transform_single_builders.hpp"
#include "HeiConnect/set_cover/solver_greedy.hpp"
#include "HeiConnect/set_cover/solver_greedy_cheapest.hpp"
#include "HeiConnect/set_cover/solver_ilp.hpp"

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
}

TEST(ConnectivityAugmentationSolvers, CSRGreedyIncreasesConnectivity)
{
    expect_csr_solver_increases_connectivity(GreedySetCoverSolver<>{});
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
