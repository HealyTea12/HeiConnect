#include <chrono>
#include <filesystem>
#include <memory>

#include "algorithm_registry.hpp"
#include "algorithm_runner.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/set_cover/common.hpp"
#include "HeiConnect/set_cover/solver_lazy_block_tree_ilp.hpp"

namespace
{
    class LazyBlockTreeILPRunner : public AlgorithmRunner
    {
    public:
        void run(
            const std::filesystem::path& graph_file,
            const std::filesystem::path& link_file,
            const std::filesystem::path&) override
        {
            auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
            auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);
            USSolution solution;
            HeiConnect::LazyBlockTreeSolverILP<> solver;

            const auto start = std::chrono::steady_clock::now();
            if (!solver.solve(graph, link_graph, solution))
            {
                return;
            }
            result.time_total = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

            double solution_cost = 0.0;
            for (const size_t link_id : solution.get_solution())
            {
                solution_cost += link_graph.weights.at(link_id);
            }
            result.solution_cost = solution_cost;
            result.solution_size = solution.get_solution_size();
        }
    };

    [[maybe_unused]] const bool registered_lazy_block_tree_ilp = [] {
        global_registry.add(
            "Lazy Block Tree ILP",
            "Exact block-tree ILP with lazy residual constraints",
            [](const ParamMap&) { return std::make_unique<LazyBlockTreeILPRunner>(); });
        return true;
    }();
}
