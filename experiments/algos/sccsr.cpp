#include <random>
#include <utility>

#include "algorithm_registry.hpp"
#include "algorithm_runner.hpp"
#include "common_stages.hpp"
#include "HeiConnect/conn_aug/reducers/full_single_dom_reducer.hpp"
#include "HeiConnect/sc_reduction/transform_single_builders.hpp"
#include "HeiConnect/set_cover/local_search/local_search.hpp"
#include "HeiConnect/set_cover/local_search/p_perturbation.hpp"
#include "HeiConnect/set_cover/solver_greedy.hpp"
#include "HeiConnect/set_cover/trimmer.hpp"
#include "HeiConnect/set_cover/util.hpp"

class SetCoverGreedyCSR : public AlgorithmRunner
{
public:
    void run(const std::filesystem::path& graph_file) override
    {
        const ReadFromFileCSRStage read_stage{};
        const GraphMetricsCalculator graph_metrics_stage{};
        FullSingleDomReducer<1> full_single_dom_reducer{true, true};
        auto build_stage = ConstructSetCoverStage{};
        auto build_context_stage = [&](auto set_cover) {
            BoundContext ctx{BasicContext(*set_cover), USSolution{}};
            return std::tuple{std::move(set_cover), std::move(ctx)};
        };
        GreedySetCoverSolver<1> greedy_solver{};
        SetCoverTrimmer<1> trimmer{};
        auto engine = std::mt19937_64{std::random_device{}()};
        auto move_generator = CostFractionPerturbator(.5, engine);
        auto local_search =
            BreakAndRepairSearch<decltype(move_generator), decltype(greedy_solver), 1, DefaultEvaluator, true>{
                move_generator,
                greedy_solver,
                5.0};

        auto pipeline = Pipeline{
            read_stage,
            graph_metrics_stage,
            full_single_dom_reducer,
            graph_metrics_stage,
            build_stage,
            build_context_stage,
            greedy_solver,
            trimmer,
            local_search};

        auto pipeline_result = pipeline.run(graph_file);
        pipeline_metrics = std::move(pipeline_result.second);

        auto final_state = std::move(pipeline_result.first);
        const auto& set_cover = *std::get<0>(final_state);
        const auto& context = std::get<1>(final_state);
        result.solution_cost = HeiConnect::sc::cost(set_cover, context.get_solution());
        result.solution_size = context.get_solution().size();
    }
};

namespace
{
    [[maybe_unused]] const bool registered_set_cover_csr = [] {
        global_registry.add("Set Cover CSR", "Set Cover algorithm using CSR representation", [](const ParamMap&) {
            return std::make_unique<SetCoverGreedyCSR>();
        });
        return true;
    }();
}
