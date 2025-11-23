#include <random>
#include <filesystem>
#include <numeric>

#include "HeiConnect/graph.hpp"
#include "HeiConnect/set_cover/transform_single.hpp"
#include "HeiConnect/greedy.hpp"
#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/tools/timer.hpp"
#include "HeiConnect/data_structures/graph_utils.hpp"
#include "dataset_manager.hpp"

void write_cycle_graphs(std::filesystem::path graph_dir)
{
    auto a = 0;
    std::vector<size_t> cycle_sizes = {};
    for (size_t n = 200; n <= 200; n += 10)
        cycle_sizes.push_back(n);
    for (size_t n_nodes : cycle_sizes)
    {
        auto cycle_graph = create_cycle_graph(n_nodes);
        std::unordered_map<size_t, std::vector<size_t>> map_to_original_graph{};
        for (size_t i = 0; i < n_nodes; ++i)
        {
            map_to_original_graph[i + 1] = {i + 1};
        }
        cycle_graph.write_to_file_graphML(graph_dir / ("cycle_" + std::to_string(n_nodes) + ".xml"), map_to_original_graph);
        cycle_graph.write_to_file_metis(graph_dir / ("cycle_" + std::to_string(n_nodes) + ".graph"));
    }
}

void experiment(std::filesystem::path graph_dir)
{
    auto curr_dir = std::filesystem::current_path();
    auto output_dir = curr_dir / "output";
    std::filesystem::create_directories(output_dir);

    for (auto file : std::filesystem::directory_iterator(graph_dir))
    {
        if (file.path().filename().string().find("cycle") == std::string::npos)
            continue;
        if (!file.path().filename().string().ends_with(".xml"))
            continue;
        std::cout << "Processing graph: " << file.path() << std::endl;
        try
        {
            WeightedCRFGraph<> graph = WeightedCRFGraph<>::read_from_file_graphML(file.path());
            if (graph.graph.vertices.size() - 1 > 200)
                continue;
            auto link_graph = graph.generate_links([](size_t u, size_t v)
                                                   { return 1.0; });

            {
                auto timer = Timer{
                    "Cycle " + std::to_string(graph.graph.vertices.size() - 1),
                    output_dir / "set_cover_greedy"};
                auto sc = construct_set_cover(
                    graph.graph.vertices,
                    graph.graph.edges,
                    graph.weights,
                    link_graph.graph.vertices,
                    link_graph.graph.edges,
                    link_graph.weights);
                timer.add_checkpoint("Reduction");
                SetCoverSolverGreedyParallel solver{sc};
                solver.solve();
                auto solution = solver.get_solution();
            }

            auto g = graph::GraphPair{};
            auto src_dir = std::filesystem::current_path().parent_path();
            g.read_graph(file.path().parent_path() / (file.path().stem().string() + ".graph"),
                         file.path());
            g.add_links(1, 1.f, 0); // all links = 1.0
            // std::cout << g.cactus.num_nodes() << " cactus nodes." << std::endl;
            // std::cout << g.cactus.num_edges() << " cactus edges." << std::endl;
            // std::cout << g.original_graph.num_nodes() << " original graph nodes." << std::endl;
            // std::cout << g.original_graph.num_edges() << " original graph edges." << std::endl;
            // std::cout << g.cactus.links[0][0].weight << " link weight." << std::endl;
            // std::cout << std::accumulate(g.cactus.links.begin(), g.cactus.links.end(), 0ull,
            //   [](size_t acc, const std::vector<graph::Edge> &vec)
            //    {
            // return acc + vec.size();
            //})
            //             << " links added." << std::endl;
            {
                auto timer = Timer{
                    "Cycle " + std::to_string(graph.graph.vertices.size() - 1),
                    output_dir / "direct_greedy"};
                auto solution = solver::greedy_heuristic_strong(g);
            };
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error processing graph " << file.path() << ": " << e.what() << std::endl;
        }
    }
}

int main()
{
    auto graph_dir = std::filesystem::temp_directory_path() / "graphs";
    std::filesystem::create_directory(graph_dir);
    write_cycle_graphs(graph_dir);
    experiment(
        graph_dir
        // std::filesystem::current_path() / "../graphs/misc"
    );
    std::filesystem::remove_all(graph_dir);
}
