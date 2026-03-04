#include <HeiConnect/set_cover/transform_single.hpp>
#include <HeiConnect/data_structures/graph_utils.hpp>
#include <benchmark/benchmark.h>

static void BM_scred_csr(benchmark::State &state)
{
    WeightedCRFGraph<> graph = create_star_graph(state.range(0));
    WeightedCRFGraph<> link_graph = graph.generate_links([](size_t u, size_t v)
                                                         { return 1.0; });
    volatile int dummy = 0;
    for (auto _ : state)
    {
        SetCover sc = construct_set_cover(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        dummy = 1;
        benchmark::ClobberMemory();
    }
}

static void BM_scred_bit(benchmark::State &state)
{
    WeightedCRFGraph<> graph = create_star_graph(state.range(0));
    WeightedCRFGraph<> link_graph = graph.generate_links([](size_t u, size_t v)
                                                         { return 1.0; });
    volatile int dummy = 0;
    for (auto _ : state)
    {
        SetCoverBit sc = construct_set_cover_bit_matrix(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        dummy = 1;
        benchmark::ClobberMemory();
    }
}

static void BM_scred_partial_bit(benchmark::State &state)
{
    WeightedCRFGraph<> graph = create_star_graph(state.range(0));
    WeightedCRFGraph<> link_graph = graph.generate_links([](size_t u, size_t v)
                                                         { return 1.0; });
    for (auto _ : state)
    {
        SetCoverPseudo sc = construct_set_cover_pseudo(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        benchmark::DoNotOptimize(sc.min_cuts.data());
        benchmark::ClobberMemory();
    }
}

static void BM_scred_partial_ancestry(benchmark::State &state)
{
    WeightedCRFGraph<> graph = create_star_graph(state.range(0));
    WeightedCRFGraph<> link_graph = graph.generate_links([](size_t u, size_t v)
                                                         { return 1.0; });
    volatile int dummy = 0;
    for (auto _ : state)
    {
        SetCoverPseudo sc = construct_set_cover_pseudo_ancestry(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        dummy = 1;
        benchmark::ClobberMemory();
    }
}

BENCHMARK(BM_scred_partial_ancestry)->DenseRange(128, 1024, 128);
BENCHMARK(BM_scred_partial_bit)->DenseRange(128, 1024, 128);
BENCHMARK(BM_scred_bit)->DenseRange(128, 1024, 128);
BENCHMARK(BM_scred_csr)->DenseRange(128, 1024, 128);
BENCHMARK_MAIN();