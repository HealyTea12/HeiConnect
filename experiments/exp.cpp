#include <random>
#include <filesystem>
#include <numeric>

#include "boost/program_options.hpp"

#include "HeiConnect/graph.hpp"
#include "HeiConnect/set_cover/transform_single.hpp"
#include "HeiConnect/greedy.hpp"
#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/tools/timer.hpp"
#include "HeiConnect/data_structures/graph_utils.hpp"
#include "HeiConnect/ilp.hpp"
#include "dataset_manager.hpp"

enum class Algorithms
{
    SetCoverGreedySingleThreadedPQ,
    SetCoverSharpGreedy,
    SetCoverILP,
    DirectGreedy,
    GWC,
    MSTConnect,
    DirectILP
};

std::array<std::string, 7> algorithm_names = {
    "SetCoverGreedySingleThreadedPQ",
    "SetCoverSharpGreedy",
    "SetCoverILP",
    "DirectGreedy",
    "GWC",
    "MSTConnect",
    "DirectILP"};

static Algorithms from_string(const std::string &algo)
{
    for (size_t i = 0; i < algorithm_names.size(); ++i)
    {
        if (algorithm_names[i] == algo)
        {
            return static_cast<Algorithms>(i);
        }
    }
    throw std::invalid_argument("Unknown algorithm: " + algo);
}

void experiment(std::filesystem::path graph_dir, std::filesystem::path output_file, Algorithms algorithm)
{

    for (auto file : std::filesystem::directory_iterator(graph_dir))
    {
        if (!file.path().filename().string().ends_with(".xml"))
            continue;
        std::cout << "Processing graph: " << file.path() << std::endl;
        try
        {
            if (algorithm == Algorithms::SetCoverGreedySingleThreadedPQ)
            {
                WeightedCRFGraph<> graph = WeightedCRFGraph<>::read_from_file_graphML(file.path());
                auto link_graph = graph.generate_links([](size_t u, size_t v)
                                                       { return 1.0; });

                {
                    auto timer = Timer{
                        file.path().filename().string(),
                        output_file};
                    auto sc = construct_set_cover(
                        graph.graph.vertices,
                        graph.graph.edges,
                        graph.weights,
                        link_graph.graph.vertices,
                        link_graph.graph.edges,
                        link_graph.weights);
                    timer.add_checkpoint("Reduction");
                    SetCoverSolverGreedySingleThreadedPQ solver{std::move(sc)};
                    solver.solve();
                    auto solution = solver.get_solution();
                    std::cout << "Solution cost: " << solver.get_solution_cost() << std::endl;
                }
            }
            else if (algorithm == Algorithms::DirectGreedy)
            {

                auto g = graph::GraphPair{};
                g.read_graph(file.path().parent_path() / (file.path().stem().string() + ".graph"),
                             file.path());
                g.add_links(1, 1.f, 0); // all links = 1.0
                {
                    auto timer = Timer{
                        file.path().filename().string(),
                        output_file};
                    auto solution = solver::greedy_heuristic_strong(g);
                    auto solution_cost = 0.0;
                    for (const auto &edge : solution)
                    {
                        solution_cost += edge.weight;
                    }
                    std::cout << "Solution cost: " << solution_cost << std::endl;
                };
            }
            else if (algorithm == Algorithms::GWC)
            {
                graph::GraphPair g;
                g.read_graph(file.path().parent_path() / (file.path().stem().string() + ".graph"),
                             file.path());
                g.add_links(1, 1.f, 0); // all links = 1.0
                graph::DynamicCactus g_dynamic;
                {
                    auto timer = Timer{
                        file.path().filename().string(),
                        output_file};
                    g_dynamic.read_from_file(file.path());
                    g_dynamic.copy_links(g);
                    auto solution = solver::greedy_dynamic_bounds(g_dynamic);
                    auto solution_cost = 0.0;
                    for (const auto &edge : solution)
                    {
                        solution_cost += edge.weight;
                    }
                    std::cout << "Solution cost: " << solution_cost << std::endl;
                }
            }
            else if (algorithm == Algorithms::DirectILP)
            {
                auto g = graph::GraphPair{};
                g.read_graph(file.path().parent_path() / (file.path().stem().string() + ".graph"),
                             file.path());
                g.add_links(1, 1.f, 0); // all links = 1.0
                {
                    auto timer = Timer{
                        file.path().filename().string(),
                        output_file};
                    auto solution = solver::ilp(g, false, 0);
                    auto solution_cost = 0.0;
                    for (const auto &edge : solution)
                    {
                        solution_cost += edge.weight;
                    }
                    std::cout << "Solution cost: " << solution_cost << std::endl;
                };
            }
            else if (algorithm == Algorithms::SetCoverILP)
            {
                WeightedCRFGraph<> graph = WeightedCRFGraph<>::read_from_file_graphML(file.path());
                auto link_graph = graph.generate_links([](size_t u, size_t v)
                                                       { return 1.0; });

                {
                    auto timer = Timer{
                        file.path().filename().string(),
                        output_file};
                    SetCover sc = construct_set_cover(
                        graph.graph.vertices,
                        graph.graph.edges,
                        graph.weights,
                        link_graph.graph.vertices,
                        link_graph.graph.edges,
                        link_graph.weights);
                    timer.add_checkpoint("Reduction");
                    // ILP solver
                    SetCoverSolverILP solver{std::move(sc)};
                    solver.solve();
                    auto ilp_solution = solver.get_solution();
                    double total_cost = 0.0;
                    for (const auto &set_index : ilp_solution)
                    {
                        total_cost += sc.costs[set_index];
                    }
                    std::cout << "Solution cost: " << total_cost << std::endl;
                }
            }
            else if (algorithm == Algorithms::MSTConnect)
            {
                auto g = graph::GraphPair{};
                g.read_graph(file.path().parent_path() / (file.path().stem().string() + ".graph"),
                             file.path());
                g.add_links(1, 1.f, 0); // all links = 1.0
                {
                    auto timer = Timer{
                        file.path().filename().string(),
                        output_file};
                    auto solution = solver::greedy_mst_max_flow(g).first;
                    auto solution_cost = 0.0;
                    for (const auto &edge : solution)
                    {
                        solution_cost += edge.weight;
                    }
                    std::cout << "Solution cost: " << solution_cost << std::endl;
                };
            }
            else if (algorithm == Algorithms::SetCoverSharpGreedy)
            {
                WeightedCRFGraph<> graph = WeightedCRFGraph<>::read_from_file_graphML(file.path());
                auto link_graph = graph.generate_links([](size_t u, size_t v)
                                                       { return 1.0; });

                {
                    auto timer = Timer{
                        file.path().filename().string(),
                        output_file};
                    SetCover sc = construct_set_cover(
                        graph.graph.vertices,
                        graph.graph.edges,
                        graph.weights,
                        link_graph.graph.vertices,
                        link_graph.graph.edges,
                        link_graph.weights);
                    timer.add_checkpoint("Reduction");
                    SetCoverSolverSharpGreedy solver{std::move(sc)};
                    solver.solve();
                    solver.trim_solution();
                    auto solution = solver.get_solution();
                    std::cout << "Solution cost: " << solver.get_solution_cost() << std::endl;
                }
            }
            std::cout << "----------------------------------------" << std::endl;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error processing graph " << file.path() << ": " << e.what() << std::endl;
        }
    }
}

namespace po = boost::program_options;

int main(int argc, char **argv)
{
    std::stringstream available_algorihtms_to_use{};
    for (const auto &name : algorithm_names)
    {
        available_algorihtms_to_use << name << ";\n";
    }
    po::options_description desc("Options");
    desc.add_options()("help", "help message")("output_file,o", po::value<std::string>())("algorithm,a", po::value<std::string>(), available_algorihtms_to_use.str().c_str())("input_dir,i", po::value<std::string>(), "input graphp directory (recursive scan for .graph and .xml files)");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help") || !vm.count("output_file") || !vm.count("algorithm") || !vm.count("input_dir"))
    {
        std::cout << desc << std::endl;
        return 1;
    }
    auto graph_dir = std::filesystem::path(vm["input_dir"].as<std::string>());
    auto output_file = std::filesystem::path(vm["output_file"].as<std::string>());
    auto algorithm_str = vm["algorithm"].as<std::string>();
    Algorithms algorithm = from_string(algorithm_str);
    experiment(graph_dir, output_file, algorithm);
}
