#include <HeiConnect/sc_reduction/transform_single_builders.hpp>
#include <HeiConnect/data_structures/graph_utils.hpp>
#include <benchmark/benchmark.h>

static void BM_scred_csr(benchmark::State &state)
{
    auto graph = create_star_graph<uint64_t, uint64_t, double>(state.range(0));
    auto link_graph = graph.generate_links([](uint64_t u, uint64_t v)
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
    auto graph = create_star_graph<uint64_t, uint64_t, double>(state.range(0));
    auto link_graph = graph.generate_links([](uint64_t u, uint64_t v)
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
    auto graph = create_star_graph<uint64_t, uint64_t, double>(state.range(0));
    auto link_graph = graph.generate_links([](uint64_t u, uint64_t v)
                                           { return 1.0; });
    for (auto _ : state)
    {
        SetCoverPseudo<uint64_t, uint64_t> sc = construct_set_cover_pseudo(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        benchmark::DoNotOptimize(sc.get_num_sets());
        benchmark::ClobberMemory();
    }
}

static void BM_scred_partial_ancestry(benchmark::State &state)
{
    auto graph = create_star_graph<uint64_t, uint64_t, double>(state.range(0));
    auto link_graph = graph.generate_links([](uint64_t u, uint64_t v)
                                           { return 1.0; });
    volatile int dummy = 0;
    for (auto _ : state)
    {
        auto sc = construct_set_cover_pseudo_ancestry(
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

static void BM_scred_partial_ancestry_vec(benchmark::State &state)
{
    auto graph = create_star_graph<uint32_t, uint32_t, double>(state.range(0));
    auto link_graph = graph.generate_links([](uint32_t u, uint32_t v)
                                           { return 1.0; });
    volatile int dummy = 0;
    for (auto _ : state)
    {
        SetCoverPseudo<uint32_t, uint32_t> sc = construct_set_cover_pseudo_ancestry_vec(
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

const int MIN = 256;
const int MAX = 1 << 10;
const int STEP = 2;

BENCHMARK(BM_scred_partial_ancestry_vec)->Range(MIN, MAX);
BENCHMARK(BM_scred_partial_ancestry)->Range(MIN, MAX);
BENCHMARK(BM_scred_partial_bit)->Range(MIN, MAX);
BENCHMARK(BM_scred_bit)->Range(MIN, MAX);
BENCHMARK(BM_scred_csr)->Range(MIN, MAX);
BENCHMARK_MAIN();