#include <chrono>
#include <cctype>
#include <random>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>

#include <omp.h>

#include "algorithm_registry.hpp"
#include "algorithm_runner.hpp"
#include "common_stages.hpp"
#include "HeiConnect/conn_aug/reducers/full_single_dom_reducer.hpp"
#include "HeiConnect/sc_reduction/transform_single_builders.hpp"
#include "HeiConnect/set_cover/local_search/local_search.hpp"
#include "HeiConnect/set_cover/local_search/p_perturbation.hpp"
#include "HeiConnect/set_cover/solver_greedy.hpp"
#include "HeiConnect/set_cover/solver_greedy_cheapest.hpp"
#include "HeiConnect/set_cover/solver_ilp.hpp"
#include "HeiConnect/set_cover/trimmer.hpp"
#include "HeiConnect/set_cover/util.hpp"
#include "HeiConnect/tools/timer.hpp"
#include "HeiConnect/visualization/dot_writer.hpp"

struct CSRReductionArtifacts
{
    std::vector<size_t> reduced_to_original;
    std::vector<double> original_costs;
    std::unordered_set<size_t> forced_original_links;
};

enum class SetCoverRepresentation
{
    CSR,
    Oracle,
    DoubleCSRCyc
};

enum class SetCoverSolver
{
    Greedy,
    GreedyCheapest,
    ILP
};

struct SetCoverGreedyCSRConfig
{
    bool run_reductions = true;
    bool draw_graphs = false;
    ConnectivityAugmentationReductionConfig reduction_config{};
    bool run_local_search = true;
    double local_search_time_seconds = 60.0;
    SetCoverRepresentation representation = SetCoverRepresentation::CSR;
    SetCoverSolver solver = SetCoverSolver::Greedy;
};

struct ConstructSetCoverOracleStage
{
    static constexpr std::string_view name = "Construct set cover (oracle)";

    template<typename GraphType, typename LinkGraphType>
    auto operator()(const GraphType& graph, const LinkGraphType& link_graph) const
    {
        auto set_cover = construct_set_cover_oracle(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        return std::make_shared<const decltype(set_cover)>(std::move(set_cover));
    }
};

struct ConstructSetCoverDoubleStage
{
    static constexpr std::string_view name = "Construct set cover (double CSR/Cyc)";

    template<typename GraphType, typename LinkGraphType>
    auto operator()(const GraphType& graph, const LinkGraphType& link_graph) const
    {
        auto set_cover = construct_set_cover_cyc_pseudo_ancestry_vec(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        return std::make_shared<const decltype(set_cover)>(std::move(set_cover));
    }
};

struct BuildBasicContextStage
{
    template<typename SetCoverType>
    auto operator()(std::shared_ptr<const SetCoverType> set_cover) const
    {
        BoundContext context{BasicContext(*set_cover), USSolution{}};
        return std::tuple{std::move(set_cover), std::move(context)};
    }
};

struct BuildDoubleContextStage
{
    template<typename SetCoverType>
    auto operator()(std::shared_ptr<const SetCoverType> set_cover) const
    {
        BitPackedContext pseudo_context{set_cover->first()};
        CycContext cyc_context{set_cover->second()};
        ContextDouble context{*set_cover, std::move(pseudo_context), std::move(cyc_context)};
        BoundContext bound_context{std::move(context), USSolution{}};
        return std::tuple{std::move(set_cover), std::move(bound_context)};
    }
};

template<typename NodeID, typename LinkEdgeID, typename LinkEdgeWeight>
WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>
extend_link_graph_vertices(const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph, size_t num_vertices)
{
    std::vector<NodeID> vertices = link_graph.graph.vertices;
    vertices.resize(num_vertices + 1, static_cast<NodeID>(link_graph.num_edges()));
    return {{std::move(vertices), link_graph.graph.edges}, link_graph.weights};
}

template<typename LocalSearchType>
class OptionalLocalSearchStage
{
public:
    OptionalLocalSearchStage(LocalSearchType& local_search, bool run_local_search) :
        m_localSearch(local_search),
        m_runLocalSearch(run_local_search)
    {}

    static constexpr std::string_view name = "Local search";

    template<typename SetCoverType, typename ContextType>
    auto operator()(std::shared_ptr<const SetCoverType> set_cover, ContextType context)
    {
        if (!m_runLocalSearch)
        {
            return std::tuple{std::move(set_cover), std::move(context)};
        }
        return m_localSearch(std::move(set_cover), std::move(context));
    }

    std::optional<StageMetrics> emit_metrics() const
    {
        if (!m_runLocalSearch)
        {
            return StageMetrics{{"skipped", "true"}};
        }
        return m_localSearch.emit_metrics();
    }

private:
    LocalSearchType& m_localSearch;
    bool m_runLocalSearch;
};

class CSRReductionStage
{
public:
    CSRReductionStage(
        std::shared_ptr<CSRReductionArtifacts> artifacts,
        bool draw_graphs,
        ConnectivityAugmentationReductionConfig reduction_config,
        std::filesystem::path before_reductions_file,
        std::filesystem::path after_reductions_file) :
        m_artifacts(std::move(artifacts)),
        m_draw_graphs(draw_graphs),
        m_reducer(reduction_config),
        m_before_reductions_file(std::move(before_reductions_file)),
        m_after_reductions_file(std::move(after_reductions_file))
    {}

    static std::string name()
    {
        return "CSR connectivity reductions";
    }

    auto operator()(const WeightedCRFGraph<>& graph, const WeightedCRFGraph<>& link_graph)
    {
        using Clock = std::chrono::high_resolution_clock;
        m_metrics = StageMetrics{};

        auto setup_start = Clock::now();
        auto link_remap = make_identity_link_remap(link_graph);
        UnionFind uf(graph.num_vertices());
        USSolution forced_solution;

        m_artifacts->original_costs = link_graph.weights;
        auto setup_end = Clock::now();
        auto before_dot_start = Clock::now();
        if (m_draw_graphs)
        {
            HeiConnect::visualization::write_to_dot(graph, link_graph, m_before_reductions_file);
        }
        auto before_dot_end = Clock::now();

        auto reducer_start = Clock::now();
        auto [reduced_graph, reduced_link_graph, final_link_remap, final_uf] =
            m_reducer.run(graph, link_graph, link_remap, uf, forced_solution);
        auto reducer_end = Clock::now();
        (void)final_uf;
        m_artifacts->forced_original_links = forced_solution.get_solution();
        auto after_dot_start = Clock::now();
        if (m_draw_graphs)
        {
            HeiConnect::visualization::write_to_dot(reduced_graph, reduced_link_graph, m_after_reductions_file);
        }
        auto after_dot_end = Clock::now();

        auto artifact_mapping_start = Clock::now();
        std::map<ConnAugLink, size_t> original_id_by_link;
        for (const auto& [link, origin] : final_link_remap)
        {
            original_id_by_link[link] = origin.original_id;
        }

        m_artifacts->reduced_to_original.resize(reduced_link_graph.num_edges());
        for (size_t u = 0; u < reduced_link_graph.num_vertices(); ++u)
        {
            for (size_t e = reduced_link_graph.graph.vertices[u]; e < reduced_link_graph.graph.vertices[u + 1]; ++e)
            {
                const size_t v = reduced_link_graph.graph.edges[e];
                m_artifacts->reduced_to_original[e] = original_id_by_link.at(normalize_link(u, v));
            }
        }
        auto artifact_mapping_end = Clock::now();

        if (auto reducer_metrics = m_reducer.emit_metrics(); reducer_metrics.has_value())
        {
            m_metrics->insert(m_metrics->end(), reducer_metrics->begin(), reducer_metrics->end());
        }
        m_metrics->push_back(
            {"reducer_total_time",
             HeiConnect::tools::format_duration(reducer_start, reducer_end, HeiConnect::tools::TimeUnit::Seconds)});
        m_metrics->push_back({"graphs_drawn", m_draw_graphs ? "true" : "false"});
        append_stage_timings(
            setup_start,
            setup_end,
            before_dot_start,
            before_dot_end,
            after_dot_start,
            after_dot_end,
            artifact_mapping_start,
            artifact_mapping_end);

        return std::tuple{std::move(reduced_graph), std::move(reduced_link_graph)};
    }

    std::optional<StageMetrics> emit_metrics() const
    {
        return m_metrics;
    }

private:
    template<typename TimePoint>
    void append_stage_timings(
        TimePoint setup_start,
        TimePoint setup_end,
        TimePoint before_dot_start,
        TimePoint before_dot_end,
        TimePoint after_dot_start,
        TimePoint after_dot_end,
        TimePoint artifact_mapping_start,
        TimePoint artifact_mapping_end)
    {
        m_metrics->push_back(
            {"stage_setup_time",
             HeiConnect::tools::format_duration(setup_start, setup_end, HeiConnect::tools::TimeUnit::Seconds)});
        m_metrics->push_back(
            {"before_reductions_dot_write_time",
             HeiConnect::tools::format_duration(
                 before_dot_start,
                 before_dot_end,
                 HeiConnect::tools::TimeUnit::Seconds)});
        m_metrics->push_back(
            {"after_reductions_dot_write_time",
             HeiConnect::tools::format_duration(after_dot_start, after_dot_end, HeiConnect::tools::TimeUnit::Seconds)});
        m_metrics->push_back(
            {"artifact_mapping_time",
             HeiConnect::tools::format_duration(
                 artifact_mapping_start,
                 artifact_mapping_end,
                 HeiConnect::tools::TimeUnit::Seconds)});
    }

    std::shared_ptr<CSRReductionArtifacts> m_artifacts;
    bool m_draw_graphs;
    std::filesystem::path m_before_reductions_file;
    std::filesystem::path m_after_reductions_file;
    FullSingleDomReducer<1> m_reducer;
    std::optional<StageMetrics> m_metrics;
};

class SetCoverGreedyCSR : public AlgorithmRunner
{
public:
    explicit SetCoverGreedyCSR(SetCoverGreedyCSRConfig config) : m_config(config)
    {}

    void print_results(std::ostream& os) override
    {
        AlgorithmRunner::print_results(os);
        if (result.solution_cost.has_value())
        {
            os << "Complete original solution cost: " << *result.solution_cost << "\n";
        }
        if (result.solution_size.has_value())
        {
            os << "Complete original solution size: " << *result.solution_size << "\n";
        }
    }

    void
    run(const std::filesystem::path& graph_file,
        const std::filesystem::path& link_file,
        const std::filesystem::path& output_dir) override
    {
        switch (m_config.representation)
        {
            case SetCoverRepresentation::CSR:
                run_with_representation(
                    ConstructSetCoverStage{},
                    BuildBasicContextStage{},
                    graph_file,
                    link_file,
                    output_dir);
                return;
            case SetCoverRepresentation::Oracle:
                run_with_representation(
                    ConstructSetCoverOracleStage{},
                    BuildBasicContextStage{},
                    graph_file,
                    link_file,
                    output_dir);
                return;
            case SetCoverRepresentation::DoubleCSRCyc:
                run_with_representation(
                    ConstructSetCoverDoubleStage{},
                    BuildDoubleContextStage{},
                    graph_file,
                    link_file,
                    output_dir);
                return;
        }
    }

private:
    template<typename BuildStage, typename BuildContextStage>
    void run_with_representation(
        BuildStage build_stage,
        BuildContextStage build_context_stage,
        const std::filesystem::path& graph_file,
        const std::filesystem::path& link_file,
        const std::filesystem::path& output_dir)
    {
        switch (m_config.solver)
        {
            case SetCoverSolver::Greedy:
                run_with_solver(
                    std::move(build_stage),
                    std::move(build_context_stage),
                    GreedySetCoverSolver<1>{},
                    graph_file,
                    link_file,
                    output_dir);
                return;
            case SetCoverSolver::GreedyCheapest:
                run_with_solver(
                    std::move(build_stage),
                    std::move(build_context_stage),
                    SetCoverSolverGreedyCheapest{},
                    graph_file,
                    link_file,
                    output_dir);
                return;
            case SetCoverSolver::ILP:
                run_with_solver(
                    std::move(build_stage),
                    std::move(build_context_stage),
                    SetCoverSolverILP<1>{},
                    graph_file,
                    link_file,
                    output_dir);
                return;
        }
    }

    template<typename BuildStage, typename BuildContextStage, typename SolverStage>
    void run_with_solver(
        BuildStage build_stage,
        BuildContextStage build_context_stage,
        SolverStage solver,
        const std::filesystem::path& graph_file,
        const std::filesystem::path& link_file,
        const std::filesystem::path& output_dir)
    {
        const ReadFromFileCSRStage read_stage{};
        const GraphMetricsCalculator graph_metrics_stage{};
        SetCoverTrimmer<1> trimmer{};
        GreedySetCoverSolver<1> repair_solver{};
        auto engine = std::mt19937_64{std::random_device{}()};
        auto move_generator = CostFractionPerturbator(.1, engine);
        auto local_search =
            BreakAndRepairSearch<decltype(move_generator), decltype(repair_solver), 1, DefaultEvaluator, true>{
                move_generator,
                repair_solver,
                m_config.local_search_time_seconds};
        const bool run_local_search = m_config.run_local_search && m_config.solver != SetCoverSolver::ILP &&
            m_config.representation != SetCoverRepresentation::DoubleCSRCyc;
        OptionalLocalSearchStage optional_local_search{local_search, run_local_search};

        if (!m_config.run_reductions)
        {
            auto pipeline = Pipeline{
                read_stage,
                graph_metrics_stage,
                build_stage,
                build_context_stage,
                solver,
                trimmer,
                optional_local_search};

            auto pipeline_result = pipeline.run(graph_file, link_file);
            pipeline_metrics = std::move(pipeline_result.second);

            auto final_state = std::move(pipeline_result.first);
            const auto& set_cover = *std::get<0>(final_state);
            const auto& context = std::get<1>(final_state);
            result.solution_cost = HeiConnect::sc::cost(set_cover, context);
            result.solution_size = context.get_solution().size();
            return;
        }

        auto reduction_artifacts = std::make_shared<CSRReductionArtifacts>();
        CSRReductionStage reduction_stage{
            reduction_artifacts,
            m_config.draw_graphs,
            m_config.reduction_config,
            output_dir / (graph_file.stem().string() + "_before_reductions.dot"),
            output_dir / (graph_file.stem().string() + "_after_reductions.dot")};
        auto pipeline = Pipeline{
            read_stage,
            graph_metrics_stage,
            reduction_stage,
            graph_metrics_stage,
            build_stage,
            build_context_stage,
            solver,
            trimmer,
            optional_local_search};

        auto pipeline_result = pipeline.run(graph_file, link_file);
        pipeline_metrics = std::move(pipeline_result.second);

        auto final_state = std::move(pipeline_result.first);
        const auto& context = std::get<1>(final_state);
        std::unordered_set<size_t> original_solution = reduction_artifacts->forced_original_links;
        for (const size_t reduced_id : context.get_solution())
        {
            original_solution.insert(reduction_artifacts->reduced_to_original.at(reduced_id));
        }

        double solution_cost = 0.0;
        for (const size_t original_id : original_solution)
        {
            solution_cost += reduction_artifacts->original_costs.at(original_id);
        }
        result.solution_cost = solution_cost;
        result.solution_size = original_solution.size();
    }

    SetCoverGreedyCSRConfig m_config;
};

class BlockTreeSetCoverGreedyCSR : public AlgorithmRunner
{
public:
    explicit BlockTreeSetCoverGreedyCSR(SetCoverGreedyCSRConfig config) : m_config(config)
    {}

    void print_results(std::ostream& os) override
    {
        AlgorithmRunner::print_results(os);
        if (result.solution_cost.has_value())
        {
            os << "Complete original solution cost: " << *result.solution_cost << "\n";
        }
        if (result.solution_size.has_value())
        {
            os << "Complete original solution size: " << *result.solution_size << "\n";
        }
    }

    void
    run(const std::filesystem::path& graph_file,
        const std::filesystem::path& link_file,
        const std::filesystem::path& output_dir) override
    {
        (void)output_dir;
        pipeline_metrics.stages.clear();

        const double read_start = omp_get_wtime();
        auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
        auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);
        record_stage(
            "Read from file (CSR)",
            read_start,
            omp_get_wtime(),
            StageMetrics{
                {"n", std::to_string(graph.num_vertices())},
                {"m", std::to_string(graph.num_edges())},
                {"num_links", std::to_string(link_graph.num_edges())},
            });
        std::unordered_set<size_t> original_solution;
        const auto original_costs = link_graph.weights;

        const double block_tree_start = omp_get_wtime();
        auto [block_tree, cycle_ids] = graph.cactus_generate_block_tree(0);
        auto block_tree_link_graph = extend_link_graph_vertices(link_graph, block_tree.num_vertices());
        record_stage(
            "Build block tree",
            block_tree_start,
            omp_get_wtime(),
            StageMetrics{
                {"n", std::to_string(block_tree.num_vertices())},
                {"m", std::to_string(block_tree.num_edges())},
                {"num_cycles", std::to_string(cycle_ids.size())},
                {"num_links", std::to_string(block_tree_link_graph.num_edges())},
            });

        std::unordered_set<size_t> first_solution;
        if (m_config.run_reductions)
        {
            const double reduction_start = omp_get_wtime();
            auto block_tree_link_remap = make_identity_link_remap(block_tree_link_graph);
            UnionFind block_tree_uf(block_tree.num_vertices());
            USSolution forced_solution;
            FullSingleDomReducer<1> reducer{m_config.reduction_config};
            auto [reduced_graph, reduced_link_graph, reduced_link_remap, ignored_uf] = reducer.run(
                block_tree,
                block_tree_link_graph,
                block_tree_link_remap,
                block_tree_uf,
                forced_solution);
            (void)ignored_uf;
            block_tree = std::move(reduced_graph);
            block_tree_link_graph = std::move(reduced_link_graph);
            block_tree_link_remap = std::move(reduced_link_remap);
            first_solution = forced_solution.get_solution();
            record_stage(
                "Block-tree reductions",
                reduction_start,
                omp_get_wtime(),
                reducer.emit_metrics().value_or(StageMetrics{}));

            const auto reduced_solution = solve_set_cover("block_tree", block_tree, block_tree_link_graph);
            for (const size_t reduced_id : reduced_solution)
            {
                const auto link = link_by_edge_id(block_tree_link_graph, reduced_id);
                first_solution.insert(block_tree_link_remap.at(link).original_id);
            }
        }
        else
        {
            first_solution = solve_set_cover("block_tree", block_tree, block_tree_link_graph);
        }
        original_solution.insert(first_solution.begin(), first_solution.end());

        const double contraction_start = omp_get_wtime();
        auto link_remap = make_identity_link_remap(link_graph);
        UnionFind contraction_uf(graph.num_vertices());
        const auto merge_stats =
            add_links_to_union_find_frozen_stack(graph, link_graph, contraction_uf, first_solution);
        auto [contracted_graph, contracted_link_graph, node_remap] =
            materialize_contractions(graph, link_graph, contraction_uf, link_remap);
        (void)node_remap;
        record_stage(
            "Contract first solution",
            contraction_start,
            omp_get_wtime(),
            StageMetrics{
                {"selected_links", std::to_string(first_solution.size())},
                {"total_merged_nodes", std::to_string(merge_stats.total_merged_nodes)},
                {"n", std::to_string(contracted_graph.num_vertices())},
                {"m", std::to_string(contracted_graph.num_edges())},
                {"num_links", std::to_string(contracted_link_graph.num_edges())},
            });

        if (contracted_graph.num_vertices() <= 1 || contracted_link_graph.num_edges() == 0)
        {
            write_original_solution_result(original_solution, original_costs);
            return;
        }

        if (m_config.run_reductions)
        {
            const double reduction_start = omp_get_wtime();
            USSolution forced_solution;
            UnionFind reduction_uf(contracted_graph.num_vertices());
            FullSingleDomReducer<1> reducer{m_config.reduction_config};
            auto [reduced_graph, reduced_link_graph, reduced_link_remap, ignored_uf] =
                reducer.run(contracted_graph, contracted_link_graph, link_remap, reduction_uf, forced_solution);
            (void)ignored_uf;
            contracted_graph = std::move(reduced_graph);
            contracted_link_graph = std::move(reduced_link_graph);
            link_remap = std::move(reduced_link_remap);
            original_solution.insert(forced_solution.get_solution().begin(), forced_solution.get_solution().end());
            record_stage(
                "Real-instance reductions",
                reduction_start,
                omp_get_wtime(),
                reducer.emit_metrics().value_or(StageMetrics{}));
        }

        if (contracted_graph.num_vertices() <= 1 || contracted_link_graph.num_edges() == 0)
        {
            write_original_solution_result(original_solution, original_costs);
            return;
        }

        const auto second_solution = solve_set_cover("contracted_real", contracted_graph, contracted_link_graph);
        for (const size_t reduced_id : second_solution)
        {
            const auto link = link_by_edge_id(contracted_link_graph, reduced_id);
            original_solution.insert(link_remap.at(link).original_id);
        }

        write_original_solution_result(original_solution, original_costs);
    }

private:
    void
    record_stage(std::string stage_name, double start, double end, std::optional<StageMetrics> metrics = std::nullopt)
    {
        PerStageMetrics stage_metrics;
        stage_metrics.duration_seconds = end - start;
        stage_metrics.metrics = std::move(metrics);
        stage_metrics.stage_name = std::move(stage_name);
        pipeline_metrics.stages.emplace_back(std::move(stage_metrics));
    }

    void write_original_solution_result(
        const std::unordered_set<size_t>& original_solution,
        const std::vector<double>& original_costs)
    {
        double solution_cost = 0.0;
        for (const size_t original_id : original_solution)
        {
            solution_cost += original_costs.at(original_id);
        }
        result.solution_cost = solution_cost;
        result.solution_size = original_solution.size();

        const double now = omp_get_wtime();
        record_stage(
            "Final original solution",
            now,
            now,
            StageMetrics{
                {"cost", std::to_string(solution_cost)},
                {"size", std::to_string(original_solution.size())},
            });
    }

    std::unordered_set<size_t> solve_set_cover(
        const std::string& prefix,
        const WeightedCRFGraph<>& graph,
        const WeightedCRFGraph<>& link_graph)
    {
        if (graph.num_vertices() <= 1 || link_graph.num_edges() == 0)
        {
            return {};
        }

        switch (m_config.representation)
        {
            case SetCoverRepresentation::CSR:
                return solve_with_representation(
                    prefix,
                    graph,
                    link_graph,
                    ConstructSetCoverStage{},
                    BuildBasicContextStage{});
            case SetCoverRepresentation::Oracle:
                return solve_with_representation(
                    prefix,
                    graph,
                    link_graph,
                    ConstructSetCoverOracleStage{},
                    BuildBasicContextStage{});
            case SetCoverRepresentation::DoubleCSRCyc:
                return solve_with_representation(
                    prefix,
                    graph,
                    link_graph,
                    ConstructSetCoverDoubleStage{},
                    BuildDoubleContextStage{});
        }
        return {};
    }

    template<typename BuildStage, typename BuildContextStage>
    std::unordered_set<size_t> solve_with_representation(
        const std::string& prefix,
        const WeightedCRFGraph<>& graph,
        const WeightedCRFGraph<>& link_graph,
        BuildStage build_stage,
        BuildContextStage build_context_stage)
    {
        const double construct_start = omp_get_wtime();
        auto set_cover = build_stage(graph, link_graph);
        record_stage(
            prefix + " " + std::string(BuildStage::name),
            construct_start,
            omp_get_wtime(),
            StageMetrics{
                {"num_sets", std::to_string(set_cover->get_num_sets())},
                {"num_elements", std::to_string(set_cover->get_num_elements())},
            });
        auto context_state = build_context_stage(std::move(set_cover));

        switch (m_config.solver)
        {
            case SetCoverSolver::Greedy:
                return solve_with_solver(prefix, std::move(context_state), GreedySetCoverSolver<1>{});
            case SetCoverSolver::GreedyCheapest:
                return solve_with_solver(prefix, std::move(context_state), SetCoverSolverGreedyCheapest{});
            case SetCoverSolver::ILP:
                return solve_with_solver(prefix, std::move(context_state), SetCoverSolverILP<1>{});
        }
        return {};
    }

    template<typename State, typename SolverStage>
    std::unordered_set<size_t> solve_with_solver(const std::string& prefix, State state, SolverStage solver)
    {
        auto set_cover = std::move(std::get<0>(state));
        auto context = std::move(std::get<1>(state));
        SetCoverTrimmer<1> trimmer{};
        GreedySetCoverSolver<1> repair_solver{};
        auto engine = std::mt19937_64{std::random_device{}()};
        auto move_generator = CostFractionPerturbator(.5, engine);
        auto local_search =
            BreakAndRepairSearch<decltype(move_generator), decltype(repair_solver), 1, DefaultEvaluator, true>{
                move_generator,
                repair_solver,
                m_config.local_search_time_seconds};

        const double solve_start = omp_get_wtime();
        solver.solve(*set_cover, context);
        std::optional<StageMetrics> solver_metrics;
        if constexpr (requires { solver.emit_metrics(); })
        {
            solver_metrics = solver.emit_metrics();
        }
        record_stage(
            prefix + " " + std::string(SolverStage::name),
            solve_start,
            omp_get_wtime(),
            std::move(solver_metrics));

        const double trim_start = omp_get_wtime();
        trimmer.trim(*set_cover, context);
        record_stage(prefix + " trim redundant sets", trim_start, omp_get_wtime(), trimmer.emit_metrics());

        const double local_search_start = omp_get_wtime();
        const bool run_local_search = m_config.run_local_search && m_config.solver != SetCoverSolver::ILP &&
            m_config.representation != SetCoverRepresentation::DoubleCSRCyc;
        if (run_local_search)
        {
            local_search.run(*set_cover, context);
            record_stage(prefix + " local search", local_search_start, omp_get_wtime(), local_search.emit_metrics());
        }
        else
        {
            record_stage(
                prefix + " local search",
                local_search_start,
                omp_get_wtime(),
                StageMetrics{{"skipped", "true"}});
        }

        return context.get_solution();
    }

    ConnAugLink link_by_edge_id(const WeightedCRFGraph<>& link_graph, size_t edge_id) const
    {
        for (size_t u = 0; u < link_graph.num_vertices(); ++u)
        {
            for (size_t e = link_graph.graph.vertices[u]; e < link_graph.graph.vertices[u + 1]; ++e)
            {
                if (e == edge_id)
                {
                    return normalize_link(u, link_graph.graph.edges[e]);
                }
            }
        }
        throw std::out_of_range("Link edge id not found");
    }

    SetCoverGreedyCSRConfig m_config;
};

namespace
{
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

    SetCoverRepresentation parse_representation_param(
        const ParamMap& params,
        SetCoverRepresentation default_value)
    {
        const auto it = params.find("reduction_type");
        if (it == params.end())
        {
            return default_value;
        }

        const std::string value = normalized_param_value(it->second);
        if (value == "csr")
        {
            return SetCoverRepresentation::CSR;
        }
        if (value == "oracle")
        {
            return SetCoverRepresentation::Oracle;
        }
        if (value == "double" || value == "double_csr_cyc" || value == "csr_cyc")
        {
            return SetCoverRepresentation::DoubleCSRCyc;
        }
        throw std::invalid_argument("Parameter 'reduction_type' must be csr, oracle, or double");
    }

    SetCoverSolver parse_solver_param(const ParamMap& params, SetCoverSolver default_value)
    {
        const auto it = params.find("solver");
        if (it == params.end())
        {
            return default_value;
        }

        const std::string value = normalized_param_value(it->second);
        if (value == "greedy")
        {
            return SetCoverSolver::Greedy;
        }
        if (value == "greedycheapest" || value == "greedy_cheapest" || value == "cheapest")
        {
            return SetCoverSolver::GreedyCheapest;
        }
        if (value == "ilp")
        {
            return SetCoverSolver::ILP;
        }
        throw std::invalid_argument("Parameter 'solver' must be greedy, greedy_cheapest, or ilp");
    }

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

    double parse_double_param(const ParamMap& params, const std::string& name, double default_value)
    {
        const auto it = params.find(name);
        if (it == params.end())
        {
            return default_value;
        }
        const double value = std::stod(it->second);
        if (value < 0.0)
        {
            throw std::invalid_argument("Parameter '" + name + "' must be non-negative");
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

    SetCoverGreedyCSRConfig parse_set_cover_greedy_csr_config(
        const ParamMap& params,
        SetCoverGreedyCSRConfig config = SetCoverGreedyCSRConfig{})
    {
        config.run_reductions = parse_bool_param(params, "reductions", config.run_reductions);
        config.draw_graphs = parse_bool_param(params, "draw_graphs", config.draw_graphs);
        config.reduction_config.run_shortest_path_reduction =
            parse_bool_param(params, "shortest_path_reduction", config.reduction_config.run_shortest_path_reduction);
        config.reduction_config.run_project_in =
            parse_bool_param(params, "project_in", config.reduction_config.run_project_in);
        config.reduction_config.run_project_out =
            parse_bool_param(params, "project_out", config.reduction_config.run_project_out);
        config.reduction_config.run_cycle_reduction =
            parse_bool_param(params, "cycle_reduction", config.reduction_config.run_cycle_reduction);
        config.reduction_config.max_rounds = parse_size_param(params, "max_rounds", config.reduction_config.max_rounds);
        config.reduction_config.compute_shortest_paths =
            parse_bool_param(params, "compute_shortest_paths", config.reduction_config.compute_shortest_paths);
        config.run_local_search = parse_bool_param(params, "local_search", config.run_local_search);
        config.run_local_search = !parse_bool_param(params, "skip_local_search", !config.run_local_search);
        config.local_search_time_seconds =
            parse_double_param(params, "local_search_time_seconds", config.local_search_time_seconds);
        config.representation = parse_representation_param(params, config.representation);
        config.solver = parse_solver_param(params, config.solver);
        return config;
    }

    [[maybe_unused]] const bool registered_set_cover_csr = [] {
        global_registry.add(
            "Set Cover CSR",
            "Configurable set-cover representation and solver",
            [](const ParamMap& params) {
                return std::make_unique<SetCoverGreedyCSR>(parse_set_cover_greedy_csr_config(params));
            });
        global_registry.add(
            "Block Tree Set Cover CSR",
            "Set Cover CSR with a block-tree contraction solve first",
            [](const ParamMap& params) {
                SetCoverGreedyCSRConfig config;
                config.run_reductions = false;
                config.run_local_search = false;
                return std::make_unique<BlockTreeSetCoverGreedyCSR>(parse_set_cover_greedy_csr_config(params, config));
            });
        return true;
    }();
}
