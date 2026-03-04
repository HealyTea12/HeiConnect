#include <benchmark/benchmark.h>
#include <vector>

static void BM_DenseRange(benchmark::State &state)
{
    for (auto _ : state)
    {
        std::vector<int> v(state.range(0), state.range(0));
        auto data = v.data();
        benchmark::DoNotOptimize(data);
        benchmark::ClobberMemory();
    }
}
BENCHMARK(BM_DenseRange)->DenseRange(0, 1024, 128);
BENCHMARK_MAIN();