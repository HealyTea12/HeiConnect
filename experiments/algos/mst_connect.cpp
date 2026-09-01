#include <chrono>
#include <filesystem>
#include <list>
#include <stdexcept>
#include <string>

#include "algorithm_registry.hpp"
#include "algorithm_runner.hpp"
#include "HeiConnect/graph.hpp"
#include "HeiConnect/greedy.hpp"

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

    int parse_nonnegative_int_param(const ParamMap& params, const std::string& name, int default_value)
    {
        const auto parameter = params.find(name);
        if (parameter == params.end())
        {
            return default_value;
        }

        size_t parsed_characters = 0;
        const int value = std::stoi(parameter->second, &parsed_characters);
        if (parsed_characters != parameter->second.size() || value < 0)
        {
            throw std::invalid_argument("Parameter '" + name + "' must be a nonnegative integer");
        }
        return value;
    }

    graph::GraphPair load_graph_pair(
        const std::filesystem::path& graph_file,
        const std::filesystem::path& link_file)
    {
        graph::GraphPair graph_pair;
        graph_pair.read_graph(graph_file.parent_path() / (graph_file.stem().string() + ".graph"), graph_file);
        graph_pair.add_links(link_file, 1.0, 0);
        return graph_pair;
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

    class MSTConnectRunner : public AlgorithmRunner
    {
    public:
        void run(
            const std::filesystem::path& graph_file,
            const std::filesystem::path& link_file,
            const std::filesystem::path&) override
        {
            auto graph_pair = load_graph_pair(graph_file, link_file);

            const auto start = std::chrono::steady_clock::now();
            const auto solution = solver::greedy_mst_max_flow(graph_pair).first;
            result.time_total = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            result.solution_cost = calculate_solution_cost(solution);
            result.solution_size = solution.size();
        }
    };

    class MSTConnectLocalSearchRunner : public AlgorithmRunner
    {
    public:
        MSTConnectLocalSearchRunner(int depth, bool cache, int trees)
            : m_depth(depth), m_cache(cache), m_trees(trees)
        {}

        void run(
            const std::filesystem::path& graph_file,
            const std::filesystem::path& link_file,
            const std::filesystem::path&) override
        {
            auto graph_pair = load_graph_pair(graph_file, link_file);

            const auto start = std::chrono::steady_clock::now();
            const auto solution =
                solver::greedy_2mst_localsearch_flow(graph_pair, m_depth, m_cache, m_trees);
            result.time_total = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            result.solution_cost = calculate_solution_cost(solution);
            result.solution_size = solution.size();
        }

    private:
        int m_depth;
        bool m_cache;
        int m_trees;
    };

    [[maybe_unused]] const bool registered_mst_connect = [] {
        global_registry.add(
            "mst-connect",
            "MST-Connect heuristic",
            [](const ParamMap&) { return std::make_unique<MSTConnectRunner>(); });

        const std::vector<AlgorithmParameter> local_search_parameters{
            {"depth", "0", "Maximum alternating-path search depth."},
            {"cache", "false", "Cache invalid alternating paths."},
            {"trees", "0", "Number of MSTs; 0 uses the algorithm default of 2."}};
        global_registry.add(
            "mst-connect-ls",
            "Local search on an MST-Connect solution",
            local_search_parameters,
            [](const ParamMap& params) {
                return std::make_unique<MSTConnectLocalSearchRunner>(
                    parse_nonnegative_int_param(params, "depth", 0),
                    parse_bool_param(params, "cache", false),
                    parse_nonnegative_int_param(params, "trees", 0));
            });
        return true;
    }();
}
