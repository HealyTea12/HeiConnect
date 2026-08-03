#include <benchmark/benchmark.h>

#include "HeiConnect/conn_aug/reducers/cycle_reducer.hpp"
#include "HeiConnect/data_structures/intersection_index/intersection_tree.hpp"

static std::vector<std::tuple<int, int, int>> fullLinkGraph(int cycle_size)
{
    std::vector<std::tuple<int, int, int>> links;
    for (int u = 0; u < cycle_size; ++u)
    {
        for (int v = u + 1; v < cycle_size; ++v)
        {
            links.emplace_back(u, v, (u + v) % cycle_size + 1);
        }
    }
    return links;
}

template<typename IntersectionIndex, bool ReuseIntersectionIndex>
static void BM_CycleReductionFullLinkGraph(benchmark::State& state)
{
    const int cycle_size = static_cast<int>(state.range(0));
    const auto links = fullLinkGraph(cycle_size);
    const IntersectionIndex intersection_index;

    for (auto _ : state)
    {
        const auto removable =
            cycle_domination_baseline(links, cycle_size, intersection_index, nullptr, ReuseIntersectionIndex);
        benchmark::DoNotOptimize(removable.data());
        benchmark::DoNotOptimize(removable.size());
    }

    state.SetComplexityN(cycle_size);
    state.counters["links"] = static_cast<double>(links.size());
}

BENCHMARK_TEMPLATE(BM_CycleReductionFullLinkGraph, BaselineIntersectionIdx<0>, true)
    ->RangeMultiplier(2)
    ->Range(8, 256)
    ->Complexity();
BENCHMARK_TEMPLATE(BM_CycleReductionFullLinkGraph, BaselineIntersectionIdx<0>, false)
    ->RangeMultiplier(2)
    ->Range(8, 256)
    ->Complexity();
BENCHMARK_TEMPLATE(BM_CycleReductionFullLinkGraph, IntersectionTreeIdx<0>, true)
    ->RangeMultiplier(2)
    ->Range(8, 256)
    ->Complexity();
BENCHMARK_TEMPLATE(BM_CycleReductionFullLinkGraph, IntersectionTreeIdx<0>, false)
    ->RangeMultiplier(2)
    ->Range(8, 256)
    ->Complexity();
