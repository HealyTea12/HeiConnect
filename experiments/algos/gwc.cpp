#include <filesystem>
#include <omp.h>
#include <stdexcept>

#include "algorithm_registry.hpp"
#include "algorithm_runner.hpp"
#include "HeiConnect/graph.hpp"
#include "HeiConnect/greedy.hpp"

namespace
{
    bool parse_sampling_param(const ParamMap& params)
    {
        const auto parameter = params.find("sampling");
        if (parameter == params.end() || parameter->second == "0")
        {
            return false;
        }
        if (parameter->second == "1")
        {
            return true;
        }
        throw std::invalid_argument("Parameter 'sampling' must be 0 or 1");
    }

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

class GWCRunner : public AlgorithmRunner
{
public:
    explicit GWCRunner(bool sampling) : m_sampling(sampling) {}

    void run(
        const std::filesystem::path& graph_file,
        const std::filesystem::path& link_file,
        const std::filesystem::path&) override
    {
        graph::GraphPair graph_pair;
        graph_pair.read_graph(graph_file.parent_path() / (graph_file.stem().string() + ".graph"), graph_file);
        graph_pair.add_links(link_file, 1.0, 0);

        graph::DynamicCactus dynamic_cactus;

        double start = omp_get_wtime();
        dynamic_cactus.read_from_file(graph_file);
        dynamic_cactus.copy_links(graph_pair);
        auto solution = m_sampling
                            ? solver::greedy_dynamic_sampling(dynamic_cactus)
                            : solver::greedy_dynamic_bounds(dynamic_cactus);
        double solving_time = omp_get_wtime() - start;

        result.solution_cost = calculate_solution_cost(solution);
        result.solution_size = solution.size();
        result.time_total = solving_time;
    }

private:
    bool m_sampling;
};

namespace
{
    const std::vector<AlgorithmParameter> gwc_parameters = {
        {"sampling", "0", "Set to 1 to enable sampled GWC; otherwise use bounds."}};

    [[maybe_unused]] const bool registered_gwc = [] {
        const auto factory = [](const ParamMap& params) {
            return std::make_unique<GWCRunner>(parse_sampling_param(params));
        };
        global_registry.add("gwc", "Greedy weight-coverage heuristic", gwc_parameters, factory);
        return true;
    }();
}
