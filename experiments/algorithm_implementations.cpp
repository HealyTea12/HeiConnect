#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include "algorithm_implementations.hpp"
#include "HeiConnect/link_reduction.hpp"
#include "HeiConnect/sc_reduction/transform_single_builders.hpp"
#include "HeiConnect/set_cover/file_writer.hpp"
#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/set_cover/trimmer.hpp"
#include "HeiConnect/set_cover/common.hpp"
#include "HeiConnect/set_cover/local_search/local_search.hpp"
#include "HeiConnect/set_cover/local_search/move_generator.hpp"
#include "HeiConnect/set_cover/local_search/exhaustive_move_generator.hpp"
#include "HeiConnect/set_cover/local_search/p_perturbation.hpp"
#include "HeiConnect/set_cover/solver_greedy_cheapest.hpp"
#include "HeiConnect/set_cover/util.hpp"
#include "HeiConnect/conn_aug/reducers/reduction.hpp"
#include "HeiConnect/conn_aug/reducers/star_reduction.hpp"
#include "HeiConnect/conn_aug/reducers/full_reducer.hpp"
#include "HeiConnect/pipeline/pipeline.hpp"
// #include "HeiConnect/set_cover/local_search/kset_move_generator.hpp"

namespace
{
    template<typename SetCoverType, typename ContextType>
    auto make_bound_context_stage()
    {
        return [](std::shared_ptr<const SetCoverType> set_cover) {
            ContextType context{*set_cover};
            BoundContext bound_context{std::move(context), USSolution{}};
            return std::tuple{std::move(set_cover), std::move(bound_context)};
        };
    }

    template<typename Solver>
    auto make_solver_stage(Solver& solver)
    {
        return [&solver](auto set_cover, auto context) {
            solver.solve(*set_cover, context);
            return std::tuple{std::move(set_cover), std::move(context)};
        };
    }
} // namespace

class GraphMetricsCalculator
{
public:
    struct StageMetrics
    {
        std::string instance;
        size_t n;
        size_t m;
        size_t d_min;
        size_t d_max;
        double d_avg;
    };
    auto run(const std::filesystem::path& graph_file)
    {
        m_metrics.instance = graph_file.string();
        auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
        m_metrics.n = graph.num_vertices();
        m_metrics.m = graph.num_edges();
        m_metrics.d_min = graph.min_degree();
        m_metrics.d_max = graph.max_degree();
        m_metrics.d_avg = graph.average_degree();
        return graph_file;
    }
    auto operator()(const std::filesystem::path& graph_file)
    {
        return run(graph_file);
    }

    StageMetrics emit_metrics() const
    {
        return m_metrics;
    }

private:
    StageMetrics m_metrics;
};

// ==================== Set Cover Greedy PQ ====================
void SetCoverGreedySingleThreadedPQRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    // StarReducer<1> star_reducer{4};
    // BasicLinkDomReducer<1> reducer{};
    FullReducer<1> full_reducer{true, true};
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

    auto pipeline = Pipeline{full_reducer, build_stage, build_context_stage, greedy_solver, trimmer, local_search};

    auto pipeline_result = pipeline.run(std::move(graph), std::move(link_graph));
    this->pipeline_metrics = std::move(pipeline_result.second);

    auto final_state = std::move(pipeline_result.first);
    auto sc = *std::get<0>(final_state);
    auto& bound_context = std::get<1>(final_state);
}

// ==================== Set Cover Greedy PQ Bit ====================
void SetCoverGreedySingleThreadedPQBitRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    using GraphType = std::decay_t<decltype(graph)>;
    using LinkGraphType = std::decay_t<decltype(link_graph)>;
    using SetCoverType = decltype(construct_set_cover_bit_matrix(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights));

    auto build_stage = [](const GraphType& graph, const LinkGraphType& link_graph) {
        auto constructed_set_cover = construct_set_cover_bit_matrix(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        return std::make_shared<const SetCoverType>(std::move(constructed_set_cover));
    };
    auto build_context_stage = make_bound_context_stage<SetCoverType, BitPackedContext<SetCoverType>>();
    GreedySetCoverSolver<1> greedy_solver{};
    SetCoverTrimmer<1> trimmer{};
    auto pipeline = Pipeline{build_stage, build_context_stage, greedy_solver, trimmer};

    auto pipeline_result = pipeline.run(std::move(graph), std::move(link_graph));
    this->pipeline_metrics = std::move(pipeline_result.second);

    auto final_state = std::move(pipeline_result.first);
    auto sc = *std::get<0>(final_state);
    auto& bound_context = std::get<1>(final_state);
    result.solution_cost = HeiConnect::sc::cost(sc, bound_context.get_solution());
    result.solution_size = bound_context.get_solution().size();
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;
    result.time_reduction = this->pipeline_metrics.stages[0].duration_seconds;
    result.time_solving = this->pipeline_metrics.stages[2].duration_seconds;
    result.time_trimming = this->pipeline_metrics.stages[3].duration_seconds;
    result.time_total = result.time_reduction.value_or(0.0) + this->pipeline_metrics.stages[1].duration_seconds +
        result.time_solving.value_or(0.0) + result.time_trimming.value_or(0.0);
}

// ==================== Set Cover Greedy PQ Pseudo ====================
void SetCoverGreedySingleThreadedPQPseudoRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    using GraphType = std::decay_t<decltype(graph)>;
    using LinkGraphType = std::decay_t<decltype(link_graph)>;
    using SetCoverType = decltype(construct_set_cover_pseudo(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights));

    auto build_stage = [](const GraphType& graph, const LinkGraphType& link_graph) {
        auto constructed_set_cover = construct_set_cover_pseudo(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        return std::make_shared<const SetCoverType>(std::move(constructed_set_cover));
    };
    auto build_context_stage = make_bound_context_stage<SetCoverType, BitPackedContext<SetCoverType>>();
    GreedySetCoverSolver<1> greedy_solver{};
    SetCoverTrimmer<1> trimmer{};
    auto pipeline = Pipeline{build_stage, build_context_stage, greedy_solver, trimmer};

    auto pipeline_result = pipeline.run(std::move(graph), std::move(link_graph));
    this->pipeline_metrics = std::move(pipeline_result.second);

    auto final_state = std::move(pipeline_result.first);
    auto sc = *std::get<0>(final_state);
    auto& bound_context = std::get<1>(final_state);
    result.solution_cost = HeiConnect::sc::cost(sc, bound_context.get_solution());
    result.solution_size = bound_context.get_solution().size();
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;
    result.time_reduction = this->pipeline_metrics.stages[0].duration_seconds;
    result.time_solving = this->pipeline_metrics.stages[2].duration_seconds;
    result.time_trimming = this->pipeline_metrics.stages[3].duration_seconds;
    result.time_total = result.time_reduction.value_or(0.0) + this->pipeline_metrics.stages[1].duration_seconds +
        result.time_solving.value_or(0.0) + result.time_trimming.value_or(0.0);
}

// ==================== Set Cover Greedy PQ Pseudo Ancestry ====================
void SCGWCPseudoAncestryRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    using GraphType = std::decay_t<decltype(graph)>;
    using LinkGraphType = std::decay_t<decltype(link_graph)>;
    using SetCoverType = decltype(construct_set_cover_pseudo_ancestry(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights));

    auto build_stage = [](const GraphType& graph, const LinkGraphType& link_graph) {
        auto constructed_set_cover = construct_set_cover_pseudo_ancestry(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        return std::make_shared<const SetCoverType>(std::move(constructed_set_cover));
    };
    auto build_context_stage = make_bound_context_stage<SetCoverType, BitPackedContext<SetCoverType>>();
    GreedySetCoverSolver<1> greedy_solver{};
    SetCoverTrimmer<1> trimmer{};
    auto pipeline = Pipeline{build_stage, build_context_stage, greedy_solver, trimmer};

    auto pipeline_result = pipeline.run(std::move(graph), std::move(link_graph));
    this->pipeline_metrics = std::move(pipeline_result.second);

    auto final_state = std::move(pipeline_result.first);
    auto sc = *std::get<0>(final_state);
    auto& bound_context = std::get<1>(final_state);
    result.solution_cost = HeiConnect::sc::cost(sc, bound_context.get_solution());
    result.solution_size = bound_context.get_solution().size();
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;
    result.time_reduction = this->pipeline_metrics.stages[0].duration_seconds;
    result.time_solving = this->pipeline_metrics.stages[2].duration_seconds;
    result.time_trimming = this->pipeline_metrics.stages[3].duration_seconds;
    result.time_total = result.time_reduction.value_or(0.0) + this->pipeline_metrics.stages[1].duration_seconds +
        result.time_solving.value_or(0.0) + result.time_trimming.value_or(0.0);
}

// ==================== Set Cover Greedy Cheapest ====================
void SetCoverGreedyCheapestRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    using SetCoverType = decltype(construct_set_cover(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights));

    auto build_stage = ConstructSetCoverStage{};
    auto build_context_stage = make_bound_context_stage<SetCoverType, BasicContext<SetCoverType>>();
    SetCoverSolverGreedyCheapest solver{};
    SetCoverTrimmer<1> trimmer{};
    auto solve_stage = make_solver_stage(solver);
    auto pipeline = Pipeline{build_stage, build_context_stage, solve_stage, trimmer};

    auto pipeline_result = pipeline.run(std::move(graph), std::move(link_graph));
    this->pipeline_metrics = std::move(pipeline_result.second);

    auto final_state = std::move(pipeline_result.first);
    auto sc = *std::get<0>(final_state);
    auto& bound_context = std::get<1>(final_state);
    result.solution_cost = HeiConnect::sc::cost(sc, bound_context.get_solution());
    result.solution_size = bound_context.get_solution().size();
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;
    result.time_reduction = this->pipeline_metrics.stages[0].duration_seconds;
    result.time_solving = this->pipeline_metrics.stages[2].duration_seconds;
    result.time_trimming = this->pipeline_metrics.stages[3].duration_seconds;
    result.time_total = result.time_reduction.value_or(0.0) + this->pipeline_metrics.stages[1].duration_seconds +
        result.time_solving.value_or(0.0) + result.time_trimming.value_or(0.0);
}

// ==================== (Set Cover Bit) Greedy Cheapest ====================
void SetCoverGreedyCheapestBitRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    using GraphType = std::decay_t<decltype(graph)>;
    using LinkGraphType = std::decay_t<decltype(link_graph)>;
    using SetCoverType = decltype(construct_set_cover_bit_matrix(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights));

    auto build_stage = [](const GraphType& graph, const LinkGraphType& link_graph) {
        auto constructed_set_cover = construct_set_cover_bit_matrix(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        return std::make_shared<const SetCoverType>(std::move(constructed_set_cover));
    };
    auto build_context_stage = make_bound_context_stage<SetCoverType, BitPackedContext<SetCoverType>>();
    SetCoverSolverGreedyCheapest solver{};
    SetCoverTrimmer<1> trimmer{};
    auto solve_stage = make_solver_stage(solver);
    auto pipeline = Pipeline{build_stage, build_context_stage, solve_stage, trimmer};

    auto pipeline_result = pipeline.run(std::move(graph), std::move(link_graph));
    this->pipeline_metrics = std::move(pipeline_result.second);

    auto final_state = std::move(pipeline_result.first);
    auto sc = *std::get<0>(final_state);
    auto& bound_context = std::get<1>(final_state);
    result.solution_cost = HeiConnect::sc::cost(sc, bound_context.get_solution());
    result.solution_size = bound_context.get_solution().size();
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;
    result.time_reduction = this->pipeline_metrics.stages[0].duration_seconds;
    result.time_solving = this->pipeline_metrics.stages[2].duration_seconds;
    result.time_trimming = this->pipeline_metrics.stages[3].duration_seconds;
    result.time_total = result.time_reduction.value_or(0.0) + this->pipeline_metrics.stages[1].duration_seconds +
        result.time_solving.value_or(0.0) + result.time_trimming.value_or(0.0);
}

// ==================== (Set Cover Pseudo) Greedy Cheapest ====================
void SetCoverPseudoGreedyCheapestRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    using GraphType = std::decay_t<decltype(graph)>;
    using LinkGraphType = std::decay_t<decltype(link_graph)>;
    using SetCoverType = decltype(construct_set_cover_pseudo(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights));

    auto build_stage = [](const GraphType& graph, const LinkGraphType& link_graph) {
        auto constructed_set_cover = construct_set_cover_pseudo(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        return std::make_shared<const SetCoverType>(std::move(constructed_set_cover));
    };
    auto build_context_stage = make_bound_context_stage<SetCoverType, BitPackedContext<SetCoverType>>();
    SetCoverSolverGreedyCheapest solver{};
    SetCoverTrimmer<1> trimmer{};
    auto solve_stage = make_solver_stage(solver);
    auto pipeline = Pipeline{build_stage, build_context_stage, solve_stage, trimmer};

    auto pipeline_result = pipeline.run(std::move(graph), std::move(link_graph));
    this->pipeline_metrics = std::move(pipeline_result.second);

    auto final_state = std::move(pipeline_result.first);
    auto sc = *std::get<0>(final_state);
    auto& bound_context = std::get<1>(final_state);
    result.solution_cost = HeiConnect::sc::cost(sc, bound_context.get_solution());
    result.solution_size = bound_context.get_solution().size();
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;
    result.time_reduction = this->pipeline_metrics.stages[0].duration_seconds;
    result.time_solving = this->pipeline_metrics.stages[2].duration_seconds;
    result.time_trimming = this->pipeline_metrics.stages[3].duration_seconds;
    result.time_total = result.time_reduction.value_or(0.0) + this->pipeline_metrics.stages[1].duration_seconds +
        result.time_solving.value_or(0.0) + result.time_trimming.value_or(0.0);
}
// ==================== Set Cover ILP ====================
void SetCoverILPRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    using GraphType = std::decay_t<decltype(graph)>;
    using LinkGraphType = std::decay_t<decltype(link_graph)>;
    using SetCoverType = decltype(construct_set_cover(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights));

    auto build_stage = ConstructSetCoverStage{};
    auto build_context_stage = make_bound_context_stage<SetCoverType, BasicContext<SetCoverType>>();
    SetCoverSolverILP solver{};
    auto solve_stage = make_solver_stage(solver);
    auto pipeline = Pipeline{build_stage, build_context_stage, solve_stage};

    auto pipeline_result = pipeline.run(std::move(graph), std::move(link_graph));
    this->pipeline_metrics = std::move(pipeline_result.second);

    auto final_state = std::move(pipeline_result.first);
    auto sc = *std::get<0>(final_state);
    auto& bound_context = std::get<1>(final_state);
    result.solution_cost = HeiConnect::sc::cost(sc, bound_context.get_solution());
    result.solution_size = bound_context.get_solution().size();
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;
    result.time_reduction = this->pipeline_metrics.stages[0].duration_seconds;
    result.time_solving = this->pipeline_metrics.stages[2].duration_seconds;
    result.time_total = result.time_reduction.value_or(0.0) + this->pipeline_metrics.stages[1].duration_seconds +
        result.time_solving.value_or(0.0);
}

// ==================== Set Cover ILP (Pseudo) ====================
void SetCoverPseudoILPRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    using GraphType = std::decay_t<decltype(graph)>;
    using LinkGraphType = std::decay_t<decltype(link_graph)>;
    using SetCoverType = decltype(construct_set_cover_pseudo_ancestry_vec(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights));

    auto build_stage = [](const GraphType& graph, const LinkGraphType& link_graph) {
        auto constructed_set_cover = construct_set_cover_pseudo_ancestry_vec(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        return std::make_shared<const SetCoverType>(std::move(constructed_set_cover));
    };
    auto build_context_stage = make_bound_context_stage<SetCoverType, BasicContext<SetCoverType>>();
    SetCoverSolverILP solver{};
    auto solve_stage = make_solver_stage(solver);
    auto pipeline = Pipeline{build_stage, build_context_stage, solve_stage};

    auto pipeline_result = pipeline.run(std::move(graph), std::move(link_graph));
    this->pipeline_metrics = std::move(pipeline_result.second);

    auto final_state = std::move(pipeline_result.first);
    auto sc = *std::get<0>(final_state);
    auto& bound_context = std::get<1>(final_state);
    result.solution_cost = HeiConnect::sc::cost(sc, bound_context.get_solution());
    result.solution_size = bound_context.get_solution().size();
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;
    result.time_reduction = this->pipeline_metrics.stages[0].duration_seconds;
    result.time_solving = this->pipeline_metrics.stages[2].duration_seconds;
    result.time_total = result.time_reduction.value_or(0.0) + this->pipeline_metrics.stages[1].duration_seconds +
        result.time_solving.value_or(0.0);
}

// ==================== Oracle Greedy Single Threaded PQ ====================
void OracleGreedySingleThreadedPQRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    using GraphType = std::decay_t<decltype(graph)>;
    using LinkGraphType = std::decay_t<decltype(link_graph)>;
    using SetCoverType = decltype(construct_set_cover_oracle(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights));

    auto build_stage = [](GraphType graph, LinkGraphType link_graph) {
        auto constructed_set_cover = construct_set_cover_oracle(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        return std::make_shared<const SetCoverType>(std::move(constructed_set_cover));
    };
    auto build_context_stage = make_bound_context_stage<SetCoverType, BasicContext<SetCoverType>>();
    GreedySetCoverSolver<1> greedy_solver{};
    auto pipeline = Pipeline{build_stage, build_context_stage, greedy_solver};

    auto pipeline_result = pipeline.run(std::move(graph), std::move(link_graph));
    this->pipeline_metrics = std::move(pipeline_result.second);

    auto final_state = std::move(pipeline_result.first);
    auto sc = *std::get<0>(final_state);
    auto& bound_context = std::get<1>(final_state);
    result.solution_cost = HeiConnect::sc::cost(sc, bound_context.get_solution());
    result.solution_size = bound_context.get_solution().size();
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;
    result.time_reduction = this->pipeline_metrics.stages[0].duration_seconds;
    result.time_solving = this->pipeline_metrics.stages[2].duration_seconds;
    result.time_total = result.time_reduction.value_or(0.0) + this->pipeline_metrics.stages[1].duration_seconds +
        result.time_solving.value_or(0.0);
}

// ==================== Cyc Greedy Single Threaded PQ ====================
void CycGreedySingleThreadedPQRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    using GraphType = std::decay_t<decltype(graph)>;
    using LinkGraphType = std::decay_t<decltype(link_graph)>;
    using SetCoverType = decltype(construct_set_cover_cyc_pseudo_ancestry_vec(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights));

    auto build_stage = [](const GraphType& graph, const LinkGraphType& link_graph) {
        auto constructed_set_cover = construct_set_cover_cyc_pseudo_ancestry_vec(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        return std::make_shared<const SetCoverType>(std::move(constructed_set_cover));
    };
    auto build_context_stage = [](std::shared_ptr<const SetCoverType> set_cover) {
        BitPackedContext pseudo_context{set_cover->first()};
        CycContext cyc_context{set_cover->second()};
        ContextDouble context{*set_cover, std::move(pseudo_context), std::move(cyc_context)};
        BoundContext bound_context{std::move(context), USSolution{}};
        return std::tuple{std::move(set_cover), std::move(bound_context)};
    };
    GreedySetCoverSolver<1> solver{};
    SetCoverTrimmer<1> trimmer{};
    auto pipeline = Pipeline{build_stage, build_context_stage, solver, trimmer};

    auto pipeline_result = pipeline.run(std::move(graph), std::move(link_graph));
    this->pipeline_metrics = std::move(pipeline_result.second);

    auto final_state = std::move(pipeline_result.first);
    auto sc = *std::get<0>(final_state);
    auto& bound_context = std::get<1>(final_state);
    result.solution_cost = HeiConnect::sc::cost(sc, bound_context.get_solution());
    result.solution_size = bound_context.get_solution().size();
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;
    result.time_reduction = this->pipeline_metrics.stages[0].duration_seconds;
    result.time_solving = this->pipeline_metrics.stages[2].duration_seconds;
    result.time_trimming = this->pipeline_metrics.stages[3].duration_seconds;
    result.time_total = result.time_reduction.value_or(0.0) + this->pipeline_metrics.stages[1].duration_seconds +
        result.time_solving.value_or(0.0) + result.time_trimming.value_or(0.0);
}

// ==================== Cyc Greedy Single Threaded PQ V2 ====================
void CycGreedySingleThreadedPQV2Runner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    using GraphType = std::decay_t<decltype(graph)>;
    using LinkGraphType = std::decay_t<decltype(link_graph)>;
    using SetCoverType = decltype(construct_set_cover_cyc_pseudo_ancestry_vec(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights));

    auto build_stage = [](const GraphType& graph, const LinkGraphType& link_graph) {
        auto constructed_set_cover = construct_set_cover_cyc_pseudo_ancestry_vec(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        return std::make_shared<const SetCoverType>(std::move(constructed_set_cover));
    };
    auto build_context_stage = [](std::shared_ptr<const SetCoverType> set_cover) {
        BitPackedContext pseudo_context{set_cover->first()};
        CycContext cyc_context{set_cover->second()};
        ContextDouble context{*set_cover, std::move(pseudo_context), std::move(cyc_context)};
        BoundContext bound_context{std::move(context), USSolution{}};
        return std::tuple{std::move(set_cover), std::move(bound_context)};
    };
    GreedySetCoverSolver<1> solver{};
    SetCoverTrimmer<1> trimmer{};
    auto pipeline = Pipeline{build_stage, build_context_stage, solver, trimmer};

    auto pipeline_result = pipeline.run(std::move(graph), std::move(link_graph));
    this->pipeline_metrics = std::move(pipeline_result.second);

    auto final_state = std::move(pipeline_result.first);
    auto sc = *std::get<0>(final_state);
    auto& bound_context = std::get<1>(final_state);
    result.solution_cost = HeiConnect::sc::cost(sc, bound_context.get_solution());
    result.solution_size = bound_context.get_solution().size();
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;
    result.time_reduction = this->pipeline_metrics.stages[0].duration_seconds;
    result.time_solving = this->pipeline_metrics.stages[2].duration_seconds;
    result.time_trimming = this->pipeline_metrics.stages[3].duration_seconds;
    result.time_total = result.time_reduction.value_or(0.0) + this->pipeline_metrics.stages[1].duration_seconds +
        result.time_solving.value_or(0.0) + result.time_trimming.value_or(0.0);
}

// ==================== Set Cover CSR Writer ====================
void SetCoverCsrWriterRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    double start = omp_get_wtime();
    auto sc_orig = construct_set_cover(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    double reduction_time = omp_get_wtime() - start;
    // discretize costs into integer bins and store as size_t in m_oracle
    auto sc_disc = sc_orig.discretize_costs<size_t>(10);
    {
        auto discretized = sc_disc;
        m_oracle = std::make_unique<SetCover<size_t, size_t, size_t>>(
            discretized.get_a(),
            discretized.get_b(),
            discretized.get_costs(),
            discretized.get_num_elements());
    }

    result.time_reduction = reduction_time;
    result.time_total = reduction_time;
}

void SetCoverCsrWriterRunner::print_results(std::ostream& os)
{
    std::cout << "Printing graph...\n";
    if (m_oracle)
    {
        SetCoverWriter::write(*m_oracle, os, SetCoverWriter::Format::DEFAULT);
        os << '\n';
    }
}

// ==================== Direct Greedy ====================
void DirectGreedyRunner::run(const std::filesystem::path& graph_file)
{
    auto g = load_graph_pair(graph_file);

    double start = omp_get_wtime();
    auto solution = solver::greedy_heuristic_strong(g);
    double solving_time = omp_get_wtime() - start;
    double solution_cost = calculate_edge_list_cost(solution);
    result.time_total = solving_time;
}

// ==================== GWC ====================
void GWCRunner::run(const std::filesystem::path& graph_file)
{
    auto g = load_graph_pair(graph_file);
    graph::DynamicCactus g_dynamic;

    double start = omp_get_wtime();
    g_dynamic.read_from_file(graph_file);
    g_dynamic.copy_links(g);
    auto solution = solver::greedy_dynamic_bounds(g_dynamic);
    double solving_time = omp_get_wtime() - start;
    double solution_cost = calculate_edge_list_cost(solution);
    result.solution_cost = solution_cost;
    result.solution_size = solution.size();
    result.time_total = solving_time;
}

// ==================== MST Connect ====================
void MSTConnectRunner::run(const std::filesystem::path& graph_file)
{
    auto g = load_graph_pair(graph_file);

    double start = omp_get_wtime();
    auto solution = solver::greedy_mst_max_flow(g).first;
    double solving_time = omp_get_wtime() - start;
    double solution_cost = calculate_edge_list_cost(solution);
    result.solution_cost = solution_cost;
    result.solution_size = solution.size();
    result.time_total = solving_time;
}

// ==================== Direct ILP ====================
void DirectILPRunner::run(const std::filesystem::path& graph_file)
{
    auto g = load_graph_pair(graph_file);

    double start = omp_get_wtime();
    auto solution = solver::ilp(g, false, 0);
    double solving_time = omp_get_wtime() - start;
    double solution_cost = calculate_edge_list_cost(solution);
    result.solution_cost = solution_cost;
    result.solution_size = solution.size();
    result.time_total = solving_time;
}

// ==================== Factory Function ====================
std::unique_ptr<AlgorithmRunner> create_algorithm_runner(Algorithms algorithm)
{
    switch (algorithm)
    {
        case Algorithms::SetCoverGreedySingleThreadedPQ:
            return std::make_unique<SetCoverGreedySingleThreadedPQRunner>();
        case Algorithms::SetCoverGreedySingleThreadedPQBit:
            return std::make_unique<SetCoverGreedySingleThreadedPQBitRunner>();
        case Algorithms::SetCoverGreedySingleThreadedPQPseudo:
            return std::make_unique<SetCoverGreedySingleThreadedPQPseudoRunner>();
        case Algorithms::SCGWCPseudoAncestry: return std::make_unique<SCGWCPseudoAncestryRunner>();
        case Algorithms::SetCoverGreedyCheapest: return std::make_unique<SetCoverGreedyCheapestRunner>();
        case Algorithms::SetCoverGreedyCheapestBit: return std::make_unique<SetCoverGreedyCheapestBitRunner>();
        case Algorithms::SetCoverPseudoGreedyCheapest: return std::make_unique<SetCoverPseudoGreedyCheapestRunner>();
        case Algorithms::SetCoverILP: return std::make_unique<SetCoverILPRunner>();
        case Algorithms::DirectGreedy: return std::make_unique<DirectGreedyRunner>();
        case Algorithms::GWC: return std::make_unique<GWCRunner>();
        case Algorithms::MSTConnect: return std::make_unique<MSTConnectRunner>();
        case Algorithms::DirectILP: return std::make_unique<DirectILPRunner>();
        case Algorithms::SetCoverPseudoILP: return std::make_unique<SetCoverPseudoILPRunner>();
        case Algorithms::OracleGreedySingleThreadedPQ: return std::make_unique<OracleGreedySingleThreadedPQRunner>();
        case Algorithms::CycGreedySingleThreadedPQ: return std::make_unique<CycGreedySingleThreadedPQRunner>();
        case Algorithms::CycGreedySingleThreadedPQV2: return std::make_unique<CycGreedySingleThreadedPQV2Runner>();
        case Algorithms::SetCoverCsrWriter: return std::make_unique<SetCoverCsrWriterRunner>();
        default: throw std::invalid_argument("Unknown algorithm");
    }
}
