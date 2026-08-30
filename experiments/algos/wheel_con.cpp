#include <chrono>
#include <filesystem>
#include <stdexcept>
#include <vector>

#include "algorithm_registry.hpp"
#include "algorithm_runner.hpp"
#include "HeiConnect/conn_aug/algorithms/wheel_con.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"

namespace
{
    class WheelConRunner : public AlgorithmRunner
    {
    public:
        void run(
            const std::filesystem::path& graph_file,
            const std::filesystem::path& link_file,
            const std::filesystem::path&) override
        {
            auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
            const auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);
            auto links = std::vector<HeiConnect::conn_aug::Link<>>{};
            links.reserve(link_graph.num_edges());
            for (const auto& [u, v, cost] : link_graph.csr_to_vec_links())
            {
                links.push_back({u, v, cost});
            }

            using Instance = HeiConnect::conn_aug::Instance<WeightedCRFGraph<>>;
            Instance instance{std::move(graph), std::move(links), graph_file.string()};
            HeiConnect::conn_aug::WheelCon solver;
            const auto start = std::chrono::steady_clock::now();
            const auto solution = solver.solve(instance);
            result.time_total = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            if (solution.status != HeiConnect::conn_aug::SolveStatus::Feasible)
            {
                throw std::runtime_error(solution.message);
            }
            result.solution_cost = solution.objective_value;
            result.solution_size = solution.selected_link_ids.size();
        }
    };

    [[maybe_unused]] const bool registered_wheel_con = [] {
        global_registry.add(
            "WheelCon",
            "Match opposite cactus leaves in Hamiltonian cycle order",
            [](const ParamMap&) { return std::make_unique<WheelConRunner>(); });
        return true;
    }();
}
