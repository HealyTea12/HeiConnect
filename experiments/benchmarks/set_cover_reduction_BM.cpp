#include <filesystem>

#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/sc_reduction/transform_single_builders.hpp"
#include "HeiConnect/tools/timer.hpp"

int main(int argc, char **argv)
{
    std::filesystem::path input_dir = std::filesystem::path(argv[1]);
    std::filesystem::path output_file = std::filesystem::path(argv[2]);
    std::string type = argv[3];
    for (const auto &file : std::filesystem::directory_iterator(input_dir))
    {
        if (!file.path().filename().string().ends_with(".xml"))
            continue;
        std::cout << "Processing graph: " << file.path() << std::endl;
        WeightedCRFGraph<> graph = WeightedCRFGraph<>::read_from_file_graphML(file.path());
        auto link_graph = graph.generate_links([](size_t u, size_t v)
                                               { return 1.0; });
        if (type == "csr")
        {
            Timer timer = Timer{
                file.path().filename().string() + " - CSR",
                output_file};
            double start = omp_get_wtime();
            auto set_cover = construct_set_cover(
                graph.graph.vertices,
                graph.graph.edges,
                graph.weights,
                link_graph.graph.vertices,
                link_graph.graph.edges,
                link_graph.weights);
            double end = omp_get_wtime();
            std::cout << "Constructing set cover in CSR form took " << (end - start) << " seconds." << std::endl;
            timer.add_checkpoint("reduction");
            GreedySetCoverSolver<SetCover<>, 0> solver;
            USSolution solution;
            solver.solve(set_cover, solution);
            timer.add_checkpoint("solving");
            double solution_cost = 0.0;
            for (const auto set_index : solution.get_solution())
            {
                solution_cost += set_cover.get_set_cost(set_index);
            }
            std::cout << "Solution cost: " << solution_cost << std::endl;
        }
        if (type == "bit")
        {
            Timer timer = Timer{
                file.path().filename().string() + " - BIT MATRIX",
                output_file};
            double start = omp_get_wtime();
            auto set_cover_bit = construct_set_cover_bit_matrix(
                graph.graph.vertices,
                graph.graph.edges,
                graph.weights,
                link_graph.graph.vertices,
                link_graph.graph.edges,
                link_graph.weights);
            double end = omp_get_wtime();
            std::cout << "Constructing set cover in BIT MATRIX form took " << (end - start) << " seconds." << std::endl;
            timer.add_checkpoint("reduction");
            GreedySetCoverSolver<SetCoverBit<>, 0> solver;
            USSolution solution;
            solver.solve(set_cover_bit, solution);
            timer.add_checkpoint("solving");
            double solution_cost = 0.0;
            for (const auto set_index : solution.get_solution())
            {
                solution_cost += set_cover_bit.get_set_cost(set_index);
            }
            std::cout << "Solution cost: " << solution_cost << std::endl;
        }
        if (type == "pseudo")
        {
            Timer timer = Timer{
                file.path().filename().string() + " - PSEUDO",
                output_file};
            double start = omp_get_wtime();
            auto set_cover_pseudo = construct_set_cover_pseudo(
                graph.graph.vertices,
                graph.graph.edges,
                graph.weights,
                link_graph.graph.vertices,
                link_graph.graph.edges,
                link_graph.weights);
            double end = omp_get_wtime();
            std::cout << "Constructing set cover in PSEUDO form took " << (end - start) << " seconds." << std::endl;
            timer.add_checkpoint("reduction");
            GreedySetCoverSolver<decltype(set_cover_pseudo), 0> solver;
            USSolution solution;
            solver.solve(set_cover_pseudo, solution);
            timer.add_checkpoint("solving");
            double solution_cost = 0.0;
            for (const auto set_index : solution.get_solution())
            {
                solution_cost += set_cover_pseudo.get_set_cost(set_index);
            }
            std::cout << "Solution cost: " << solution_cost << std::endl;
        }
        if (type == "pseudo_ancestry")
        {
            Timer timer = Timer{
                file.path().filename().string() + " - PSEUDO ANCESTRY",
                output_file};
            double start = omp_get_wtime();
            auto set_cover_pseudo = construct_set_cover_pseudo_ancestry(
                graph.graph.vertices,
                graph.graph.edges,
                graph.weights,
                link_graph.graph.vertices,
                link_graph.graph.edges,
                link_graph.weights);
            double end = omp_get_wtime();
            std::cout << "Constructing set cover in PSEUDO ANCESTRY form took " << (end - start) << " seconds." << std::endl;
            timer.add_checkpoint("reduction");
            GreedySetCoverSolver<decltype(set_cover_pseudo), 0> solver;
            USSolution solution;
            solver.solve(set_cover_pseudo, solution);
            timer.add_checkpoint("solving");
            double solution_cost = 0.0;
            for (const auto set_index : solution.get_solution())
            {
                solution_cost += set_cover_pseudo.get_set_cost(set_index);
            }
            std::cout << "Solution cost: " << solution_cost << std::endl;
        }
    }
}