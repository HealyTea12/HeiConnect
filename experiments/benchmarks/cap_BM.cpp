// Connecticity augmentation benchmark

#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/set_cover/transform_single.hpp"
#include <benchmark/benchmark.h>
#include <HeiConnect/data_structures/graph_utils.hpp>

static void BM_cap_pam_stars(benchmark::State &state)
{
    // Load graph and links from file
    auto g = create_star_graph<uint64_t, uint64_t, double>(state.range(0));
    WeightedCRFGraph<> link_g = g.generate_links([](size_t u, size_t v)
                                                 { return 1.0; });
    for (auto _ : state)
    {
        auto set_cover = construct_set_cover_pseudo_ancestry(
            g.graph.vertices,
            g.graph.edges,
            g.weights,
            link_g.graph.vertices,
            link_g.graph.edges,
            link_g.weights);
        SetCoverSolverGreedySingleThreadedPQ<decltype(set_cover)> solver{set_cover};
        solver.solve();
    }
}

static void BM_cap_oracle_stars(benchmark::State &state)
{
    auto g = create_star_graph<uint64_t, uint64_t, double>(state.range(0));
    WeightedCRFGraph<> link_g = g.generate_links([](size_t u, size_t v)
                                                 { return 1.0; });
    for (auto _ : state)
    {
        auto set_cover = construct_set_cover_oracle<>(
            g.graph.vertices,
            g.graph.edges,
            g.weights,
            link_g.graph.vertices,
            link_g.graph.edges,
            link_g.weights);
        SetCoverSolverGreedySingleThreadedPQ<decltype(set_cover)> solver{set_cover};
        solver.solve();
    }
}

BENCHMARK(BM_cap_pam_stars)->Range(1 << 8, 1 << 12);
BENCHMARK(BM_cap_oracle_stars)->Range(1 << 8, 1 << 12);