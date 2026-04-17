#include <benchmark/benchmark.h>

#include <filesystem>
#include <string>

#include "HeiConnect/data_structures/graph_utils.hpp"
#include "HeiConnect/greedy.hpp"
#include "HeiConnect/graph.hpp"
#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/set_cover/transform_single_builders.hpp"

static void BM_setcovercyc_cover_count(benchmark::State &state)
{
    const auto n_nodes = static_cast<int>(state.range(0));
    const auto xml_file = std::filesystem::path("datasets/cycles") / ("cycle_" + std::to_string(n_nodes) + ".xml");
    const auto links_file = xml_file.parent_path() / (xml_file.stem().string() + ".links");

    auto graph = WeightedCRFGraph<>::read_from_file_graphML(xml_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(links_file);

    auto sc = construct_set_cover_cyc_pseudo_ancestry_vec(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);

    SetCoverSolverGreedySingleThreadedPQ<decltype(sc)> solver{std::move(sc)};

    const size_t n_sets = link_graph.weights.size();
    const size_t set_index = 0;

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(solver.cover_count(set_index));
    }

    state.SetItemsProcessed(state.iterations());
    state.counters["n_sets"] = static_cast<double>(n_sets);
}

static void BM_setcovercyc_cover_count_all_sets(benchmark::State &state)
{
    const auto n_nodes = static_cast<int>(state.range(0));
    const auto xml_file = std::filesystem::path("datasets/cycles") / ("cycle_" + std::to_string(n_nodes) + ".xml");
    const auto links_file = xml_file.parent_path() / (xml_file.stem().string() + ".links");

    auto graph = WeightedCRFGraph<>::read_from_file_graphML(xml_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(links_file);

    auto sc = construct_set_cover_cyc_pseudo_ancestry_vec(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);

    SetCoverSolverGreedySingleThreadedPQ<decltype(sc)> solver{std::move(sc)};

    const size_t n_sets = link_graph.weights.size();

    for (auto _ : state)
    {
        size_t total = 0;
        for (size_t set_index = 0; set_index < n_sets; ++set_index)
        {
            total += solver.cover_count(set_index);
        }
        benchmark::DoNotOptimize(total);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations() * n_sets);
    state.counters["n_sets"] = static_cast<double>(n_sets);
}

static void BM_gwc_greedy_dynamic_bounds(benchmark::State &state)
{
    const auto n_nodes = static_cast<int>(state.range(0));
    const auto xml_file = std::filesystem::path("datasets/cycles") / ("cycle_" + std::to_string(n_nodes) + ".xml");
    const auto graph_file = xml_file.parent_path() / (xml_file.stem().string() + ".graph");
    const auto links_file = xml_file.parent_path() / (xml_file.stem().string() + ".links");

    graph::GraphPair graph_pair;
    graph_pair.read_graph(graph_file, xml_file);
    graph_pair.add_links(links_file, 1.0, 0);

    graph::DynamicCactus prototype;
    prototype.read_from_file(xml_file);
    prototype.copy_links(graph_pair);

    for (auto _ : state)
    {
        state.PauseTiming();
        graph::DynamicCactus instance = prototype;
        state.ResumeTiming();

        auto solution = solver::greedy_dynamic_bounds(instance);
        benchmark::DoNotOptimize(solution);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
}

BENCHMARK(BM_setcovercyc_cover_count)->DenseRange(100, 400, 100);
BENCHMARK(BM_setcovercyc_cover_count_all_sets)->DenseRange(100, 400, 100);
BENCHMARK(BM_gwc_greedy_dynamic_bounds)->DenseRange(100, 400, 100);
