#include <fstream>
#include <iostream>
#include <numeric>
#include "algorithm_implementations.hpp"
#include "HeiConnect/link_reduction.hpp"
#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/set_cover/trimmer.hpp"

// ==================== Set Cover Greedy PQ ====================
void SetCoverGreedySingleThreadedPQRunner::run(const std::filesystem::path &graph_file)
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
    GreedySetCoverPipeline<decltype(sc)> solver{};
    solver.solve(std::move(sc));
    double solving_time = omp_get_wtime() - start - reduction_time;

    result.solution_cost = solver.get_solution_cost();
    result.solution_size = solver.get_solution().size();
    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
}

// ==================== Set Cover Greedy PQ Bit ====================
void SetCoverGreedySingleThreadedPQBitRunner::run(const std::filesystem::path &graph_file)
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
    GreedySetCoverSolver<SetCoverBit, USSolution, BitPackedContext<SetCoverBit, USSolution>> greedy_solver{};
    BitPackedContext<SetCoverBit, USSolution> context{std::make_shared<const SetCoverBit>(sc)};
    SetCoverTrimmer<SetCoverBit, USSolution, BitPackedTrimmerContext<SetCoverBit, USSolution>> trimmer{};
    BitPackedTrimmerContext<SetCoverBit, USSolution> trimmer_context{std::make_shared<const SetCoverBit>(sc)};
    greedy_solver.solve(sc, solution, context);
    trimmer.trim(sc, solution, trimmer_context);

    // SetCoverPipeline<
    //     SetCoverBit,
    //     USSolution,
    //     SolverStageAdapter<
    //         GreedySetCoverSolver<
    //             SetCoverBit, USSolution, BitPackedContext<SetCoverBit, USSolution>>,
    //         Fresh<BitPackedContext<SetCoverBit, USSolution>>,
    //         BitPackedContext<SetCoverBit, USSolution>>,
    //     TrimmerStageAdapter<SetCoverTrimmer<SetCoverBit, USSolution, BitPackedContext<SetCoverBit, USSolution>>, Reuse, BitPackedContext<SetCoverBit, USSolution>>>
    //     solver{};
    // solver.solve(std::move(sc));
    double solving_time = omp_get_wtime() - start - reduction_time;

    result.solution_cost = std::accumulate(solution.get_solution().begin(), solution.get_solution().end(), 0.0, [&sc](double acc, size_t set_index)
                                           { return acc + sc.get_set_cost(set_index); });
    result.solution_size = solution.get_solution().size();
    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
}

// ==================== Set Cover Greedy PQ Pseudo ====================
void SetCoverGreedySingleThreadedPQPseudoRunner::run(const std::filesystem::path &graph_file)
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
    GreedySetCoverSolver<SetCoverPseudo<size_t, size_t>, USSolution, BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>> greedy_solver{};
    BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution> context{std::make_shared<const SetCoverPseudo<size_t, size_t>>(sc)};
    SetCoverTrimmer<SetCoverPseudo<size_t, size_t>, USSolution, BitPackedTrimmerContext<SetCoverPseudo<size_t, size_t>, USSolution>> trimmer{};
    BitPackedTrimmerContext<SetCoverPseudo<size_t, size_t>, USSolution> trimmer_context{std::make_shared<const SetCoverPseudo<size_t, size_t>>(sc)};
    greedy_solver.solve(sc, solution, context);
    trimmer.trim(sc, solution, trimmer_context);
    // SetCoverPipeline<
    //     SetCoverPseudo<size_t, size_t>,
    //     USSolution,
    //     SolverStageAdapter<
    //         GreedySetCoverSolver<
    //             SetCoverPseudo<size_t, size_t>, USSolution, BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>>,
    //         Fresh<BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>>,
    //         BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>>,
    //     TrimmerStageAdapter<SetCoverTrimmer<SetCoverPseudo<size_t, size_t>, USSolution, BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>>, Reuse, BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>>>
    //     solver{};
    // solver.solve(std::move(sc));
    double solving_time = omp_get_wtime() - start - reduction_time;

    result.solution_cost = std::accumulate(solution.get_solution().begin(), solution.get_solution().end(), 0.0, [&sc](double acc, size_t set_index)
                                           { return acc + sc.get_set_cost(set_index); });
    result.solution_size = solution.get_solution().size();
    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
    result.time_trimming = 0.0;
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;
    result.time_ls = 0.0;
    result.solution_cost_ls = result.solution_cost;
    result.solution_size_ls = result.solution_size;
}

// ==================== Set Cover Greedy PQ Pseudo Ancestry ====================
void SCGWCPseudoAncestryRunner::run(const std::filesystem::path &graph_file)
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
    GreedySetCoverSolver<SetCoverPseudo<size_t, size_t>, USSolution, BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>> greedy_solver{};
    BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution> context{std::make_shared<const SetCoverPseudo<size_t, size_t>>(sc)};
    SetCoverTrimmer<SetCoverPseudo<size_t, size_t>, USSolution, BitPackedTrimmerContext<SetCoverPseudo<size_t, size_t>, USSolution>> trimmer{};
    BitPackedTrimmerContext<SetCoverPseudo<size_t, size_t>, USSolution> trimmer_context{std::make_shared<const SetCoverPseudo<size_t, size_t>>(sc)};
    greedy_solver.solve(sc, solution, context);
    trimmer.trim(sc, solution, trimmer_context);
    // SetCoverPipeline<
    //     SetCoverPseudo<size_t, size_t>,
    //     USSolution,
    //     SolverStageAdapter<
    //         GreedySetCoverSolver<
    //             SetCoverPseudo<size_t, size_t>, USSolution, BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>>,
    //         Fresh<BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>>,
    //         BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>>,
    //     TrimmerStageAdapter<SetCoverTrimmer<SetCoverPseudo<size_t, size_t>, USSolution, BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>>, Reuse, BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>>>
    //     solver{};
    double solving_time = omp_get_wtime() - start - reduction_time;

    result.solution_cost = std::accumulate(solution.get_solution().begin(), solution.get_solution().end(), 0.0, [&sc](double acc, size_t set_index)
                                           { return acc + sc.get_set_cost(set_index); });
    result.solution_size = solution.get_solution().size();
    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
    result.time_trimming = 0.0;
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;
    result.time_ls = 0.0;
    result.solution_cost_ls = result.solution_cost;
    result.solution_size_ls = result.solution_size;
}

// ==================== Set Cover Greedy Cheapest ====================
void SetCoverGreedyCheapestRunner::run(const std::filesystem::path &graph_file)
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
    // CheapestSetCoverPipeline<SetCover> solver{};
    USSolution solution{};
    BasicContext<SetCover, USSolution> context{std::make_shared<const SetCover>(sc)};
    SetCoverSolverGreedyCheapest<SetCover, BasicContext<SetCover, USSolution>, USSolution> solver{};
    solver.solve(sc, solution, context);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = std::accumulate(solution.get_solution().begin(), solution.get_solution().end(), 0.0, [&sc](double acc, size_t set_index)
                                           { return acc + sc.get_set_cost(set_index); });
    result.solution_size = solution.get_solution().size();
    SetCoverTrimmer<SetCover, USSolution, BasicContext<SetCover, USSolution>> trimmer{};
    trimmer.trim(sc, solution, context);
    double trimming_time = omp_get_wtime() - start - reduction_time - solving_time;

    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time + trimming_time;

    result.solution_cost_trimmed = std::accumulate(solution.get_solution().begin(), solution.get_solution().end(), 0.0, [&sc](double acc, size_t set_index)
                                                   { return acc + sc.get_set_cost(set_index); });
    result.solution_size_trimmed = solution.get_solution().size();
}

// ==================== (Set Cover Bit) Greedy Cheapest ====================
void SetCoverGreedyCheapestBitRunner::run(const std::filesystem::path &graph_file)
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
    BitPackedContext<SetCoverBit, USSolution> context{std::make_shared<const SetCoverBit>(sc)};
    SetCoverSolverGreedyCheapest<SetCoverBit, BitPackedContext<SetCoverBit, USSolution>, USSolution> solver{};
    solver.solve(sc, solution, context);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = std::accumulate(solution.get_solution().begin(), solution.get_solution().end(), 0.0, [&sc](double acc, size_t set_index)
                                           { return acc + sc.get_set_cost(set_index); });
    result.solution_size = solution.get_solution().size();
    SetCoverTrimmer<SetCoverBit, USSolution, BitPackedTrimmerContext<SetCoverBit, USSolution>> trimmer{};
    BitPackedTrimmerContext<SetCoverBit, USSolution> trimmer_context{std::make_shared<const SetCoverBit>(sc)};
    trimmer.trim(sc, solution, trimmer_context);
    double trimming_time = omp_get_wtime() - start - reduction_time - solving_time;

    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time + trimming_time;
    result.solution_cost_trimmed = std::accumulate(solution.get_solution().begin(), solution.get_solution().end(), 0.0, [&sc](double acc, size_t set_index)
                                                   { return acc + sc.get_set_cost(set_index); });
    result.solution_size_trimmed = solution.get_solution().size();
}

// ==================== (Set Cover Pseudo) Greedy Cheapest ====================
void SetCoverPseudoGreedyCheapestRunner::run(const std::filesystem::path &graph_file)
{
    double start = omp_get_wtime();
    const std::string link_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".links");
    // const std::string original_graph_file = graph_file.parent_path() / (graph_file.filename().stem().string() + ".graph");
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
    // auto mapping = read_mapping(graph_file, read_nodes_in_original_graph(original_graph_file));
    auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);
    // auto link_reduction_output = generate_link_mapping(
    //    link_graph,
    //    graph,
    //    mapping);
    // link_graph = link_reduction_output.cactus_links;
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
    BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution> context{std::make_shared<const SetCoverPseudo<size_t, size_t>>(sc)};
    SetCoverSolverGreedyCheapest<SetCoverPseudo<size_t, size_t>, BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>, USSolution> solver{};
    solver.solve(sc, solution, context);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = std::accumulate(solution.get_solution().begin(), solution.get_solution().end(), 0.0, [&sc](double acc, size_t set_index)
                                           { return acc + sc.get_set_cost(set_index); });
    result.solution_size = solution.get_solution().size();
    SetCoverTrimmer<SetCoverPseudo<size_t, size_t>, USSolution, BitPackedTrimmerContext<SetCoverPseudo<size_t, size_t>, USSolution>> trimmer{};
    BitPackedTrimmerContext<SetCoverPseudo<size_t, size_t>, USSolution> trimmer_context{std::make_shared<const SetCoverPseudo<size_t, size_t>>(sc)};
    trimmer.trim(sc, solution, trimmer_context);
    double trimming_time = omp_get_wtime() - start - reduction_time - solving_time;

    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time + trimming_time;
    result.solution_cost_trimmed = std::accumulate(solution.get_solution().begin(), solution.get_solution().end(), 0.0, [&sc](double acc, size_t set_index)
                                                   { return acc + sc.get_set_cost(set_index); });
    result.solution_size_trimmed = solution.get_solution().size();
}
// ==================== Set Cover ILP ====================
void SetCoverILPRunner::run(const std::filesystem::path &graph_file)
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
    SetCoverSolverILP<SetCover, USSolution> solver{};
    solver.solve(sc, solution);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = std::accumulate(solution.get_solution().begin(), solution.get_solution().end(), 0.0, [&sc](double acc, size_t set_index)
                                           { return acc + sc.get_set_cost(set_index); });
    result.solution_size = solution.get_solution().size();

    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
}

// ==================== Set Cover ILP (Pseudo) ====================
void SetCoverPseudoILPRunner::run(const std::filesystem::path &graph_file)
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
    SetCoverSolverILP<SetCoverPseudo<size_t, size_t>, USSolution> solver{};
    solver.solve(sc, solution);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = std::accumulate(solution.get_solution().begin(), solution.get_solution().end(), 0.0, [&sc](double acc, size_t set_index)
                                           { return acc + sc.get_set_cost(set_index); });
    result.solution_size = solution.get_solution().size();
    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
}

// ==================== Oracle Greedy Single Threaded PQ ====================
void OracleGreedySingleThreadedPQRunner::run(const std::filesystem::path &graph_file)
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
    GreedySetCoverPipeline<decltype(sc)> solver{};
    solver.solve(std::move(sc));
    double solving_time = omp_get_wtime() - start - reduction_time;

    result.solution_cost = solver.get_solution_cost();
    result.solution_size = solver.get_solution().size();
    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
}

// ==================== Cyc Greedy Single Threaded PQ ====================
void CycGreedySingleThreadedPQRunner::run(const std::filesystem::path &graph_file)
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
    ContextDouble<
        SetCoverPseudo<size_t, size_t>,
        SetCoverCyc<size_t, size_t, double>,
        BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>,
        BasicContext<SetCoverCyc<size_t, size_t, double>, USSolution>,
        USSolution>
        context{std::make_shared<decltype(sc)>(sc)};
    GreedySetCoverSolver<decltype(sc), USSolution, decltype(context)> solver{};
    solver.solve(sc, solution, context);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = std::accumulate(solution.get_solution().begin(), solution.get_solution().end(), 0.0, [&sc](double acc, size_t set_index)
                                           { return acc + sc.get_set_cost(set_index); });
    result.solution_size = solution.get_solution().size();
    result.solution_cost_trimmed = result.solution_cost;
    result.solution_size_trimmed = result.solution_size;

    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
}

// ==================== Cyc Greedy Single Threaded PQ V2 ====================
void CycGreedySingleThreadedPQV2Runner::run(const std::filesystem::path &graph_file)
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
    ContextDouble<
        SetCoverPseudo<size_t, size_t>,
        SetCoverCyc<size_t, size_t, double>,
        BitPackedContext<SetCoverPseudo<size_t, size_t>, USSolution>,
        CycContext<size_t, size_t, double, USSolution>,
        USSolution>
        context{std::make_shared<decltype(sc)>(sc)};
    GreedySetCoverSolver<decltype(sc), USSolution, decltype(context)> solver{};
    solver.solve(sc, solution, context);
    double solving_time = omp_get_wtime() - start - reduction_time;
    result.solution_cost = std::accumulate(solution.get_solution().begin(), solution.get_solution().end(), 0.0, [&sc](double acc, size_t set_index)
                                           { return acc + sc.get_set_cost(set_index); });
    result.solution_size = solution.get_solution().size();
    // Trimmer doesn't work with this context

    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;
}

// ==================== Direct Greedy ====================
void DirectGreedyRunner::run(const std::filesystem::path &graph_file)
{
    auto g = load_graph_pair(graph_file);

    double start = omp_get_wtime();
    auto solution = solver::greedy_heuristic_strong(g);
    double solving_time = omp_get_wtime() - start;
    double solution_cost = calculate_edge_list_cost(solution);
    result.solution_cost = solution_cost;
    result.solution_size = solution.size();
    result.time_total = solving_time;
}

// ==================== GWC ====================
void GWCRunner::run(const std::filesystem::path &graph_file)
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
void MSTConnectRunner::run(const std::filesystem::path &graph_file)
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
void DirectILPRunner::run(const std::filesystem::path &graph_file)
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
    case Algorithms::SCGWCPseudoAncestry:
        return std::make_unique<SCGWCPseudoAncestryRunner>();
    case Algorithms::SetCoverGreedyCheapest:
        return std::make_unique<SetCoverGreedyCheapestRunner>();
    case Algorithms::SetCoverGreedyCheapestBit:
        return std::make_unique<SetCoverGreedyCheapestBitRunner>();
    case Algorithms::SetCoverPseudoGreedyCheapest:
        return std::make_unique<SetCoverPseudoGreedyCheapestRunner>();
    case Algorithms::SetCoverILP:
        return std::make_unique<SetCoverILPRunner>();
    case Algorithms::DirectGreedy:
        return std::make_unique<DirectGreedyRunner>();
    case Algorithms::GWC:
        return std::make_unique<GWCRunner>();
    case Algorithms::MSTConnect:
        return std::make_unique<MSTConnectRunner>();
    case Algorithms::DirectILP:
        return std::make_unique<DirectILPRunner>();
    case Algorithms::SetCoverPseudoILP:
        return std::make_unique<SetCoverPseudoILPRunner>();
    case Algorithms::OracleGreedySingleThreadedPQ:
        return std::make_unique<OracleGreedySingleThreadedPQRunner>();
    case Algorithms::CycGreedySingleThreadedPQ:
        return std::make_unique<CycGreedySingleThreadedPQRunner>();
    case Algorithms::CycGreedySingleThreadedPQV2:
        return std::make_unique<CycGreedySingleThreadedPQV2Runner>();
    default:
        throw std::invalid_argument("Unknown algorithm");
    }
}
