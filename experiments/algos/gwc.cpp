#include <filesystem>
#include <omp.h>

#include "algorithm_registry.hpp"
#include "algorithm_runner.hpp"
#include "HeiConnect/graph.hpp"
#include "HeiConnect/greedy.hpp"

namespace
{
    double calculate_solution_cost(const std::list<graph::Edge>& solution)
    {
        double cost = 0.0;
        for (const auto& edge : solution)
        {
            cost += edge.weight;
        }
        return cost;
    }
} // namespace

class ClassicalGWCRunner : public AlgorithmRunner
{
public:
    void run(const std::filesystem::path& graph_file) override
    {
        graph::GraphPair graph_pair;
        graph_pair.read_graph(graph_file.parent_path() / (graph_file.stem().string() + ".graph"), graph_file);
        graph_pair.add_links(graph_file.parent_path() / (graph_file.stem().string() + ".links"), 1.0, 0);

        graph::DynamicCactus dynamic_cactus;

        double start = omp_get_wtime();
        dynamic_cactus.read_from_file(graph_file);
        dynamic_cactus.copy_links(graph_pair);
        auto solution = solver::greedy_dynamic_bounds(dynamic_cactus);
        double solving_time = omp_get_wtime() - start;

        result.solution_cost = calculate_solution_cost(solution);
        result.solution_size = solution.size();
        result.time_total = solving_time;
    }
};

namespace
{
    [[maybe_unused]] const bool registered_classical_gwc = [] {
        global_registry.add("GWC", "Classical GWC algorithm", [](const ParamMap&) {
            return std::make_unique<ClassicalGWCRunner>();
        });
        return true;
    }();
}
