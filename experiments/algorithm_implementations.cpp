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
#include "HeiConnect/conn_aug/reduction.hpp"
#include "HeiConnect/pipeline/pipeline.hpp"
// #include "HeiConnect/set_cover/local_search/kset_move_generator.hpp"


// ==================== Set Cover Greedy PQ ====================
void SetCoverGreedySingleThreadedPQRunner::run(const std::filesystem::path& graph_file)
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
    using BoundCtx = BoundContext<BasicContext<SetCoverType>, USSolution>;

    BasicLinkDomReducer<1> reducer{};
    auto build_stage = ConstructSetCoverStage{};
    GreedySetCoverSolver<1> greedy_solver{};
    SetCoverTrimmer<1> trimmer{};
    auto engine = std::mt19937_64{std::random_device{}()};
    auto move_generator = CostFractionPerturbator(.5, engine);
    auto local_search =
        BreakAndRepairSearch<decltype(move_generator), decltype(greedy_solver), 1, DefaultEvaluator, true>{
            move_generator,
            greedy_solver,
            5.0};

    auto build_context_stage = [&](std::shared_ptr<const SetCoverType> set_cover) {
        BoundCtx ctx{BasicContext<SetCoverType>(*set_cover), USSolution{}};
        return std::tuple{std::move(set_cover), std::move(ctx)};
    };
    auto pipeline = Pipeline{reducer, build_stage, build_context_stage, greedy_solver, trimmer, local_search};

    auto pipeline_result = pipeline.run(std::move(graph), std::move(link_graph));
    this->pipeline_metrics = std::move(pipeline_result.second);

    auto final_state = std::move(pipeline_result.first);
    auto sc = *std::get<0>(final_state);
    std::cout << "Solving...\n";
    auto& bound_context = std::get<1>(final_state);
    result.solution_cost = HeiConnect::sc::cost(sc, bound_context.get_solution());
    result.solution_size = bound_context.get_solution().size();
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;
    result.solution_cost_ls = result.solution_cost;
    result.solution_size_ls = result.solution_size;
    std::cout << "Solved!\n";
}

// ==================== Set Cover Greedy PQ Bit ====================
void SetCoverGreedySingleThreadedPQBitRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    double start = omp_get_wtime();
    auto sc = construct_set_cover_bit_matrix(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    double reduction_time = omp_get_wtime() - start;
    USSolution solution{};
    BitPackedContext context{sc};
    BoundContext bound_context{context, solution};
    GreedySetCoverSolver<1> greedy_solver{};
    SetCoverTrimmer<1> trimmer{};
    greedy_solver.solve(sc, bound_context);
    result.solution_cost = HeiConnect::sc::cost(sc, solution);
    result.solution_size = solution.get_solution().size();
    trimmer.trim(sc, bound_context);
    double solving_time = omp_get_wtime() - start - reduction_time;

    result.solution_cost_trimmed = HeiConnect::sc::cost(sc, solution);
    result.solution_size_trimmed = solution.get_solution().size();
    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
}

// ==================== Set Cover Greedy PQ Pseudo ====================
void SetCoverGreedySingleThreadedPQPseudoRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    double start = omp_get_wtime();
    auto sc = construct_set_cover_pseudo(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    double reduction_time = omp_get_wtime() - start;
    USSolution solution{};
    BitPackedContext context{sc};
    BoundContext bound_context{context, solution};
    GreedySetCoverSolver<1> greedy_solver{};
    SetCoverTrimmer<1> trimmer{};
    greedy_solver.solve(sc, bound_context);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = HeiConnect::sc::cost(sc, solution);
    result.solution_size = solution.get_solution().size();
    trimmer.trim(sc, bound_context);
    double trimming_time = omp_get_wtime() - start - reduction_time - solving_time;

    result.solution_cost_trimmed = HeiConnect::sc::cost(sc, solution);
    result.solution_size_trimmed = solution.get_solution().size();
    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_trimming = trimming_time;
    result.time_total = reduction_time + solving_time + trimming_time;
}

// ==================== Set Cover Greedy PQ Pseudo Ancestry ====================
void SCGWCPseudoAncestryRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    double start = omp_get_wtime();
    auto sc = construct_set_cover_pseudo_ancestry(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    double reduction_time = omp_get_wtime() - start;
    USSolution solution{};
    BitPackedContext context{sc};
    BoundContext bound_context{context, solution};
    GreedySetCoverSolver<1> greedy_solver{};
    SetCoverTrimmer<1> trimmer{};
    greedy_solver.solve(sc, bound_context);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = HeiConnect::sc::cost(sc, solution);
    result.solution_size = solution.get_solution().size();
    trimmer.trim(sc, bound_context);
    double trimming_time = omp_get_wtime() - start - reduction_time - solving_time;

    result.solution_cost_trimmed = HeiConnect::sc::cost(sc, solution);
    result.solution_size_trimmed = solution.get_solution().size();
    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_trimming = trimming_time;
    result.time_total = reduction_time + solving_time + trimming_time;
}

// ==================== Set Cover Greedy Cheapest ====================
void SetCoverGreedyCheapestRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    double start = omp_get_wtime();
    auto sc = construct_set_cover(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    double reduction_time = omp_get_wtime() - start;
    USSolution solution{};
    BasicContext context{sc};
    BoundContext bound_context{context, solution};
    SetCoverSolverGreedyCheapest solver{};
    SetCoverTrimmer<1> trimmer{};
    solver.solve(sc, bound_context);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = HeiConnect::sc::cost(sc, solution);
    result.solution_size = solution.get_solution().size();
    trimmer.trim(sc, bound_context);
    double trimming_time = omp_get_wtime() - start - reduction_time - solving_time;

    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_trimming = trimming_time;
    result.time_total = reduction_time + solving_time + trimming_time;

    result.solution_cost_trimmed = HeiConnect::sc::cost(sc, solution);
    result.solution_size_trimmed = solution.get_solution().size();
}

// ==================== (Set Cover Bit) Greedy Cheapest ====================
void SetCoverGreedyCheapestBitRunner::run(const std::filesystem::path& graph_file)
{
    double start = omp_get_wtime();
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);
    double end = omp_get_wtime();
    std::cout << "Initialization time: " << end - start << "s\n";

    start = omp_get_wtime();
    auto sc = construct_set_cover_bit_matrix(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    double reduction_time = omp_get_wtime() - start;
    USSolution solution{};
    BitPackedContext context{sc};
    std::cout << "pruning\n";
    BoundContext bound_context{context, solution};
    SetCoverSolverGreedyCheapest solver{};
    SetCoverTrimmer<1> trimmer{};
    solver.solve(sc, bound_context);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = HeiConnect::sc::cost(sc, solution);
    result.solution_size = solution.get_solution().size();
    trimmer.trim(sc, bound_context);
    double trimming_time = omp_get_wtime() - start - reduction_time - solving_time;

    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_trimming = trimming_time;
    result.time_total = reduction_time + solving_time + trimming_time;
    result.solution_cost_trimmed = HeiConnect::sc::cost(sc, solution);
    result.solution_size_trimmed = solution.get_solution().size();
}

// ==================== (Set Cover Pseudo) Greedy Cheapest ====================
void SetCoverPseudoGreedyCheapestRunner::run(const std::filesystem::path& graph_file)
{
    double start = omp_get_wtime();
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);
    double end = omp_get_wtime();
    std::cout << "Initialization time: " << end - start << "s\n";

    start = omp_get_wtime();
    auto sc = construct_set_cover_pseudo(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    double reduction_time = omp_get_wtime() - start;
    USSolution solution{};
    BitPackedContext context{sc};
    BoundContext bound_context{context, solution};
    SetCoverSolverGreedyCheapest solver{};
    SetCoverTrimmer<1> trimmer{};
    solver.solve(sc, bound_context);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = HeiConnect::sc::cost(sc, solution);
    result.solution_size = solution.get_solution().size();
    trimmer.trim(sc, bound_context);
    double trimming_time = omp_get_wtime() - start - reduction_time - solving_time;

    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_trimming = trimming_time;
    result.time_total = reduction_time + solving_time + trimming_time;
    result.solution_cost_trimmed = HeiConnect::sc::cost(sc, solution);
    result.solution_size_trimmed = solution.get_solution().size();
}
// ==================== Set Cover ILP ====================
void SetCoverILPRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    double start = omp_get_wtime();
    SetCover sc = construct_set_cover(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    double reduction_time = omp_get_wtime() - start;
    USSolution solution{};
    BasicContext context{sc};
    BoundContext bound_context{context, solution};
    SetCoverSolverILP solver{};
    solver.solve(sc, bound_context);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = HeiConnect::sc::cost(sc, solution);
    result.solution_size = solution.get_solution().size();

    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
}

// ==================== Set Cover ILP (Pseudo) ====================
void SetCoverPseudoILPRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    double start = omp_get_wtime();
    auto sc = construct_set_cover_pseudo_ancestry_vec(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    double reduction_time = omp_get_wtime() - start;
    USSolution solution{};
    BasicContext context{sc};
    BoundContext bound_context{context, solution};
    SetCoverSolverILP solver{};
    solver.solve(sc, bound_context);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = HeiConnect::sc::cost(sc, solution);
    result.solution_size = solution.get_solution().size();
    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
}

// ==================== Oracle Greedy Single Threaded PQ ====================
void OracleGreedySingleThreadedPQRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    double start = omp_get_wtime();
    auto sc = construct_set_cover_oracle(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    double reduction_time = omp_get_wtime() - start;
    USSolution solution{};
    BasicContext context{sc};
    BoundContext bound_context{context, solution};
    GreedySetCoverSolver<1> greedy_solver{};
    greedy_solver.solve(sc, bound_context);
    double solving_time = omp_get_wtime() - start - reduction_time;

    result.solution_cost = HeiConnect::sc::cost(sc, solution);
    result.solution_size = solution.get_solution().size();
    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
}

// ==================== Cyc Greedy Single Threaded PQ ====================
void CycGreedySingleThreadedPQRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    double start = omp_get_wtime();
    auto sc = construct_set_cover_cyc_pseudo_ancestry_vec(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    double reduction_time = omp_get_wtime() - start;
    USSolution solution{};
    BitPackedContext pseudo_context{sc.first()};
    CycContext cyc_context{sc.second()};
    ContextDouble context{sc, pseudo_context, cyc_context};
    BoundContext<decltype(context), decltype(solution)> bound_context{context, solution};
    GreedySetCoverSolver<1> solver{};
    SetCoverTrimmer<1> trimmer{};
    solver.solve(sc, bound_context);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = HeiConnect::sc::cost(sc, solution);
    result.solution_size = solution.get_solution().size();
    trimmer.trim(sc, bound_context);
    result.solution_cost_trimmed = HeiConnect::sc::cost(sc, solution);
    result.solution_size_trimmed = solution.get_solution().size();

    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
}

// ==================== Cyc Greedy Single Threaded PQ V2 ====================
void CycGreedySingleThreadedPQV2Runner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    double start = omp_get_wtime();
    auto sc = construct_set_cover_cyc_pseudo_ancestry_vec(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    double reduction_time = omp_get_wtime() - start;

    USSolution solution{};
    BitPackedContext pseudo_context{sc.first()};
    CycContext cyc_context{sc.second()};
    ContextDouble context{sc, pseudo_context, cyc_context};
    BoundContext<decltype(context), decltype(solution)> bound_context{context, solution};
    GreedySetCoverSolver<1> solver{};
    SetCoverTrimmer<1> trimmer{};
    solver.solve(sc, bound_context);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = HeiConnect::sc::cost(sc, solution);
    result.solution_size = solution.get_solution().size();

    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
}

// ==================== Set Cover CSR Writer ====================
void SetCoverCsrWriterRunner::run(const std::filesystem::path& graph_file)
{
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);

    double start = omp_get_wtime();
    m_oracle = std::make_unique<SetCoverOracle<size_t, size_t, double>>(construct_set_cover_oracle(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights));
    double reduction_time = omp_get_wtime() - start;

    result.time_reduction = reduction_time;
    result.time_total = reduction_time;
}

void SetCoverCsrWriterRunner::print_results(std::ostream& os)
{
    if (m_oracle)
    {
        SetCoverWriter::write(*m_oracle, os);
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
