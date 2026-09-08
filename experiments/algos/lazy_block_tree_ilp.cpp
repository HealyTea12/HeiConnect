#include <cctype>
#include <chrono>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>

#include "algorithm_registry.hpp"
#include "algorithm_runner.hpp"
#include "HeiConnect/conn_aug/reducers/full_single_dom_reducer.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/set_cover/common.hpp"
#include "HeiConnect/set_cover/solver_lazy_block_tree_ilp.hpp"
#include "HeiConnect/visualization/dot_writer.hpp"

namespace
{
    bool parse_bool_param(const ParamMap& params, const std::string& name, bool default_value)
    {
        const auto it = params.find(name);
        if (it == params.end())
        {
            return default_value;
        }
        if (it->second == "true" || it->second == "1")
        {
            return true;
        }
        if (it->second == "false" || it->second == "0")
        {
            return false;
        }
        throw std::invalid_argument("Parameter '" + name + "' must be true or false");
    }

    std::string normalized_param_value(std::string value)
    {
        for (char& c : value)
        {
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            if (c == '-' || c == '/' || c == ' ')
            {
                c = '_';
            }
        }
        return value;
    }

    size_t parse_size_param(const ParamMap& params, const std::string& name, size_t default_value)
    {
        const auto it = params.find(name);
        if (it == params.end())
        {
            return default_value;
        }
        return std::stoull(it->second);
    }

    ConnectivityAugmentationReductionConfig parse_reduction_config(const ParamMap& params, bool& draw_graphs)
    {
        ConnectivityAugmentationReductionConfig config{};
        draw_graphs = parse_bool_param(params, "draw_graphs", false);
        config.run_project_in = parse_bool_param(params, "project_in", config.run_project_in);
        config.run_project_out = parse_bool_param(params, "project_out", config.run_project_out);
        config.run_cycle_reduction = parse_bool_param(params, "cycle_reduction", config.run_cycle_reduction);
        config.run_single_link = parse_bool_param(params, "single_link", config.run_single_link);
        config.run_element_domination = parse_bool_param(params, "element_domination", config.run_element_domination);
        config.max_rounds = parse_size_param(params, "max_rounds", config.max_rounds);
        config.compute_shortest_paths = parse_bool_param(params, "compute_shortest_paths", config.compute_shortest_paths);

        const auto it = params.find("intersection_index");
        if (it == params.end())
        {
            return config;
        }

        const std::string value = normalized_param_value(it->second);
        if (value == "baseline")
        {
            config.intersection_index = IntersectionIndexType::BASELINE;
            return config;
        }
        if (value == "intersection_tree")
        {
            config.intersection_index = IntersectionIndexType::INTERSECTION_TREE;
            return config;
        }
        if (value == "weighted_intersection_tree")
        {
            config.intersection_index = IntersectionIndexType::WEIGHTED_INTERSECTION_TREE;
            return config;
        }
        throw std::invalid_argument("Parameter 'intersection_index' must be baseline, intersection_tree, or weighted_intersection_tree");
    }

    std::unordered_set<size_t> map_solution_to_original(
        const USSolution& solution,
        const WeightedCRFGraph<>& reduced_link_graph,
        const ConnAugLinkRemap& link_remap)
    {
        std::unordered_set<size_t> original_solution;
        const auto links = reduced_link_graph.csr_to_vec_links();
        for (const size_t link_id : solution.get_solution())
        {
            const auto [u, v, _] = links.at(link_id);
            const auto it = link_remap.find(normalize_link(u, v));
            if (it != link_remap.end())
            {
                original_solution.insert(it->second.original_id);
            }
        }
        return original_solution;
    }

    double compute_solution_cost(const std::vector<double>& link_weights, const std::unordered_set<size_t>& solution)
    {
        double cost = 0.0;
        for (const size_t link_id : solution)
        {
            cost += link_weights.at(link_id);
        }
        return cost;
    }

    class LazyBlockTreeILPRunner : public AlgorithmRunner
    {
    public:
        LazyBlockTreeILPRunner(
            bool run_reductions,
            bool draw_graphs,
            ConnectivityAugmentationReductionConfig reduction_config) :
            m_run_reductions(run_reductions),
            m_draw_graphs(draw_graphs),
            m_reduction_config(std::move(reduction_config))
        {}

        void run(
            const std::filesystem::path& graph_file,
            const std::filesystem::path& link_file,
            const std::filesystem::path& output_dir) override
        {
            const auto original_graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
            const auto original_link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

            auto graph = original_graph;
            auto link_graph = original_link_graph;
            ConnAugLinkRemap link_remap;
            USSolution forced_solution;
            if (m_run_reductions)
            {
                FullSingleDomReducer<0> reducer(m_reduction_config);
                UnionFind solution_uf(original_graph.num_vertices());
                if (m_draw_graphs)
                {
                    HeiConnect::visualization::write_to_dot(graph, link_graph, output_dir / "lazy_block_tree_ilp_before_reductions.dot");
                }
                std::tie(graph, link_graph, link_remap, solution_uf) =
                    reducer.run(graph, link_graph, link_remap, solution_uf, forced_solution);
                if (m_draw_graphs)
                {
                    HeiConnect::visualization::write_to_dot(graph, link_graph, output_dir / "lazy_block_tree_ilp_after_reductions.dot");
                }
                (void)solution_uf;
            }

            HeiConnect::LazyBlockTreeSolverILP<> solver;
            USSolution solution;
            const auto start = std::chrono::steady_clock::now();
            const bool solved = solver.solve(graph, link_graph, solution);
            result.solver_status = SetCoverSolverILP<>::grb_get_status_string(solver.get_status());
            result.solution_optimal = solver.is_optimal();
            if (!solved)
            {
                return;
            }
            result.time_total = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

            if (m_run_reductions)
            {
                auto mapped_solution = map_solution_to_original(solution, link_graph, link_remap);
                for (const size_t link_id : forced_solution.get_solution())
                {
                    mapped_solution.insert(link_id);
                }

                result.solution_cost = compute_solution_cost(original_link_graph.weights, mapped_solution);
                result.solution_size = static_cast<size_t>(mapped_solution.size());
            }
    else
            {
                double cost = 0.0;
                for (const size_t link_id : solution.get_solution())
                {
                    cost += link_graph.weights.at(link_id);
                }
                result.solution_cost = cost;
                result.solution_size = solution.get_solution_size();
            }
        }

    private:
        bool m_run_reductions;
        bool m_draw_graphs;
        ConnectivityAugmentationReductionConfig m_reduction_config;
    };

    const std::vector<AlgorithmParameter> lazy_block_tree_ilp_parameters = {
        {"reductions", "true", "Enable or disable all reduction phases before solving."},
        {"draw_graphs", "false", "Write pre-reduction and post-reduction graphs as DOT files when true."},
        {"project_in", "true", "Enable projection-in reduction."},
        {"project_out", "true", "Enable projection-out reduction."},
        {"cycle_reduction", "true", "Enable cycle-reduction preprocessing."},
        {"single_link", "true", "Enable single-link reduction."},
        {"element_domination", "true", "Enable element domination reduction."},
        {"max_rounds", "0", "Maximum number of reduction rounds (0 means unlimited)."},
        {"compute_shortest_paths", "false", "Enable shortest-path computations used by intersection heuristics."},
        {"intersection_index", "baseline", "Intersection index used by cycle reduction: baseline, intersection_tree, weighted_intersection_tree."},
    };

    [[maybe_unused]] const bool registered_lazy_block_tree_ilp = [] {
        global_registry.add(
            "Lazy Block Tree ILP",
            "Exact block-tree ILP with lazy residual constraints",
            lazy_block_tree_ilp_parameters,
            [](const ParamMap& params) {
                bool draw_graphs = false;
                const auto reduction_config = parse_reduction_config(params, draw_graphs);
                return std::make_unique<LazyBlockTreeILPRunner>(
                    parse_bool_param(params, "reductions", true),
                    draw_graphs,
                    reduction_config);
            });
        return true;
    }();
}
