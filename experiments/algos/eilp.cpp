#include <chrono>
#include <filesystem>
#include <list>
#include <stdexcept>
#include <string>

#include "algorithm_registry.hpp"
#include "algorithm_runner.hpp"
#include "HeiConnect/graph.hpp"
#include "HeiConnect/ilp.hpp"

namespace
{
    bool parse_bool_param(const ParamMap& params, const std::string& name, bool default_value)
    {
        const auto parameter = params.find(name);
        if (parameter == params.end())
        {
            return default_value;
        }
        if (parameter->second == "true" || parameter->second == "1")
        {
            return true;
        }
        if (parameter->second == "false" || parameter->second == "0")
        {
            return false;
        }
        throw std::invalid_argument("Parameter '" + name + "' must be true or false");
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

    class EILPRunner : public AlgorithmRunner
    {
    public:
        explicit EILPRunner(bool use_initial) : m_use_initial(use_initial) {}

        void run(
            const std::filesystem::path& graph_file,
            const std::filesystem::path& link_file,
            const std::filesystem::path&) override
        {
            graph::GraphPair graph_pair;
            const auto original_graph = graph_file.parent_path() / (graph_file.stem().string() + ".graph");
            graph_pair.read_graph(original_graph, graph_file);
            graph_pair.add_links(link_file, 1.0, 0);

            const auto start = std::chrono::steady_clock::now();
            const auto solution = solver::ilp(graph_pair, m_use_initial);
            result.time_total = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            result.solution_cost = calculate_solution_cost(solution);
            result.solution_size = solution.size();
        }

    private:
        bool m_use_initial;
    };

    [[maybe_unused]] const bool registered_eilp = [] {
        const std::vector<AlgorithmParameter> parameters{
            {"use-initial", "false", "Use mst-connect as the initial ILP solution."}};
        const auto factory = [](const ParamMap& params) {
            return std::make_unique<EILPRunner>(parse_bool_param(params, "use-initial", false));
        };
        global_registry.add("eilp", "Optimal ILP", parameters, factory);
        return true;
    }();
}
