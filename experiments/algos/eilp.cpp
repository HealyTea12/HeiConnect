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

    int parse_int_param(const ParamMap& params, const std::string& name, int default_value)
    {
        const auto parameter = params.find(name);
        if (parameter == params.end())
        {
            return default_value;
        }

        size_t parsed_characters = 0;
        const int value = std::stoi(parameter->second, &parsed_characters);
        if (parsed_characters != parameter->second.size())
        {
            throw std::invalid_argument("Parameter '" + name + "' must be an integer");
        }
        return value;
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

    class OriginalEILPRunner : public AlgorithmRunner
    {
    public:
        OriginalEILPRunner(bool use_initial, int presolve) : m_use_initial(use_initial), m_presolve(presolve) {}

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
            const auto solution = solver::ilp(graph_pair, m_use_initial, m_presolve);
            result.time_total = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            result.solution_cost = calculate_solution_cost(solution);
            result.solution_size = solution.size();
        }

    private:
        bool m_use_initial;
        int m_presolve;
    };

    [[maybe_unused]] const bool registered_original_eilp = [] {
        global_registry.add(
            "Original EILP",
            "Original EILP solver used by the legacy app",
            std::vector<AlgorithmParameter>{
                {"use_initial", "false", "Start ILP from an initial solution when true."},
                {"presolve", "0", "Presolve level passed to the ILP solver."}},
            [](const ParamMap& params) {
                return std::make_unique<OriginalEILPRunner>(
                    parse_bool_param(params, "use_initial", false),
                    parse_int_param(params, "presolve", 0));
            });
        return true;
    }();
}
