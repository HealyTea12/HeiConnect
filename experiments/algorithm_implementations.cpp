#include <fstream>
#include <iostream>
#include "algorithm_implementations.hpp"

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
    SetCoverSolverGreedySingleThreadedPQ<SetCover> solver{std::move(sc)};
    solver.solve();
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
    SetCoverSolverGreedySingleThreadedPQ<SetCoverBit> solver{std::move(sc)};
    solver.solve();
    double solving_time = omp_get_wtime() - start - reduction_time;

    result.solution_cost = solver.get_solution_cost();
    result.solution_size = solver.get_solution().size();
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
    SetCoverSolverGreedySingleThreadedPQ<SetCoverPseudo> solver{std::move(sc)};
    solver.solve();
    double solving_time = omp_get_wtime() - start - reduction_time;

    result.solution_cost = solver.get_solution_cost();
    result.solution_size = solver.get_solution().size();
    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;

    solver.trim_solution();
    double trim_time = omp_get_wtime();
    result.time_trimming = trim_time - (start + reduction_time + solving_time);
    result.solution_cost_trimmed = solver.get_solution_cost();
    result.solution_size_trimmed = solver.get_solution().size();

    /*
    while (solver.local_search(2))
    {
        double ls_improve_time = omp_get_wtime();
        std::cout << "Improved solution: " << solver.get_solution_cost() << ',' << solver.get_solution().size() << "\n";
        std::cout << "Time so far: " << ls_improve_time - trim_time << "s\n";
        trim_time = ls_improve_time;
    }
    */
    double ls_time = omp_get_wtime();
    result.time_ls = ls_time - trim_time;
    result.solution_cost_ls = solver.get_solution_cost();
    result.solution_size_ls = solver.get_solution().size();
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
    SetCoverSolverGreedyCheapest<SetCover> solver{std::move(sc)};
    solver.solve();
    double solving_time = omp_get_wtime() - start - reduction_time;

    result.solution_cost = solver.get_solution_cost();
    result.solution_size = solver.get_solution().size();
    result.time_reduction = reduction_time;
    result.time_solving = solving_time;
    result.time_total = reduction_time + solving_time;

    solver.trim_solution();
    double trim_time = omp_get_wtime();
    result.time_trimming = trim_time - (start + reduction_time + solving_time);
    result.solution_cost_trimmed = solver.get_solution_cost();
    result.solution_size_trimmed = solver.get_solution().size();
}

// ==================== (Set Cover Pseudo) Greedy Cheapest ====================
void SetCoverPseudoGreedyCheapestRunner::run(const std::filesystem::path &graph_file)
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
    SetCoverSolverGreedyCheapest<SetCoverPseudo> solver{std::move(sc)};
    solver.solve();
    double solving_time = omp_get_wtime() - start - reduction_time;

    result.solution_cost = solver.get_solution_cost();
    result.solution_size = solver.get_solution().size();
    result.time_reduction = reduction_time;
    result.time_solving = solving_time;

    solver.trim_solution();
    double trim_time = omp_get_wtime();
    result.time_trimming = trim_time - (start + reduction_time + solving_time);
    result.solution_cost_trimmed = solver.get_solution_cost();
    result.solution_size_trimmed = solver.get_solution().size();
    result.time_total = reduction_time + solving_time + result.time_trimming;
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
    SetCoverSolverILP<SetCover> solver{std::move(sc)};
    solver.solve();
    double solving_time = omp_get_wtime() - start - reduction_time;
    auto ilp_solution = solver.get_solution();
    result.solution_cost = solver.get_solution_cost();
    result.solution_size = ilp_solution.size();
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
    case Algorithms::SetCoverGreedyCheapest:
        return std::make_unique<SetCoverGreedyCheapestRunner>();
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
    default:
        throw std::invalid_argument("Unknown algorithm");
    }
};
