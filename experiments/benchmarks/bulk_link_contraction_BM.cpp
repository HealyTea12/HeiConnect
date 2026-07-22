#include <benchmark/benchmark.h>

#include <numeric>
#include <tuple>
#include <vector>

#include "HeiConnect/conn_aug/reducers/common.hpp"
#include "HeiConnect/data_structures/graph_utils.hpp"

struct BulkContractionInstance
{
    WeightedCRFGraph<> graph;
    WeightedCRFGraph<> link_graph;
    std::vector<size_t> selected_links;
    std::vector<size_t> parent;
    std::vector<size_t> depth;
    std::vector<std::vector<int>> cycle_positions;
};

BulkContractionInstance make_crossing_cycle_instance(size_t cycle_size)
{
    auto graph = create_cycle_graph_undirected(cycle_size);
    std::vector<std::tuple<size_t, size_t, double>> links;
    links.reserve(cycle_size / 2);
    for (size_t u = 0; u < cycle_size / 2; ++u)
    {
        links.emplace_back(u, u + cycle_size / 2, 1.0);
    }

    auto link_graph = WeightedCRFGraph<>::vec_links_to_csr(links, cycle_size);
    std::vector<size_t> selected_links(links.size());
    std::iota(selected_links.begin(), selected_links.end(), 0);
    auto [block_tree, cycle_positions] = graph.cactus_generate_block_tree(0);
    auto [parent, depth] = block_tree.graph.rooted_parent_depth();
    return {
        std::move(graph),
        std::move(link_graph),
        std::move(selected_links),
        std::move(parent),
        std::move(depth),
        std::move(cycle_positions)};
}

bool is_fully_contracted(const BulkContractionInstance& instance, bool frozen_stack)
{
    UnionFind uf(instance.graph.num_vertices());
    if (frozen_stack)
    {
        add_links_to_union_find_frozen_stack(
            instance.graph.num_vertices(),
            instance.link_graph.graph,
            uf,
            instance.selected_links,
            instance.parent,
            instance.depth,
            instance.cycle_positions);
    }
    else
    {
        add_links_to_union_find(
            instance.graph.num_vertices(),
            instance.link_graph.graph,
            uf,
            instance.selected_links,
            instance.parent,
            instance.depth,
            instance.cycle_positions);
    }

    for (size_t node = 1; node < instance.graph.num_vertices(); ++node)
    {
        if (uf.find(0) != uf.find(node))
        {
            return false;
        }
    }
    return true;
}

static void BM_BulkLinkContractionCurrent(benchmark::State& state)
{
    const auto instance = make_crossing_cycle_instance(static_cast<size_t>(state.range(0)));
    if (!is_fully_contracted(instance, false))
    {
        state.SkipWithError("current contraction produced an unexpected partition");
        return;
    }
    for (auto _ : state)
    {
        UnionFind uf(instance.graph.num_vertices());
        auto stats = add_links_to_union_find(
            instance.graph.num_vertices(),
            instance.link_graph.graph,
            uf,
            instance.selected_links,
            instance.parent,
            instance.depth,
            instance.cycle_positions);
        benchmark::DoNotOptimize(stats.total_merged_nodes);
    }
    state.SetComplexityN(state.range(0));
}

static void BM_BulkLinkContractionFrozenStack(benchmark::State& state)
{
    const auto instance = make_crossing_cycle_instance(static_cast<size_t>(state.range(0)));
    if (!is_fully_contracted(instance, true))
    {
        state.SkipWithError("frozen-stack contraction produced an unexpected partition");
        return;
    }
    for (auto _ : state)
    {
        UnionFind uf(instance.graph.num_vertices());
        auto stats = add_links_to_union_find_frozen_stack(
            instance.graph.num_vertices(),
            instance.link_graph.graph,
            uf,
            instance.selected_links,
            instance.parent,
            instance.depth,
            instance.cycle_positions);
        benchmark::DoNotOptimize(stats.total_merged_nodes);
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(BM_BulkLinkContractionCurrent)->RangeMultiplier(2)->Range(64, 4096)->Complexity();
BENCHMARK(BM_BulkLinkContractionFrozenStack)->RangeMultiplier(2)->Range(64, 4096)->Complexity();
