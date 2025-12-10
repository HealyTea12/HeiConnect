#include <filesystem>
#include <gtest/gtest.h>

#include "HeiConnect/greedy.hpp"
#include "HeiConnect/data_structures/graph_utils.hpp"

TEST(GraphAugmentation, GreedyGWC)
{
    auto cycle_graph = create_cycle_graph_undirected(5); // 0 -> 1 -> 2 -> 3 ... -> 0
    std::unordered_map<size_t, std::vector<size_t>> map_to_original_graph{};
    for (size_t i{0}; i < 5; i++)
    {
        map_to_original_graph[i] = {i};
    }
    auto tmp = std::filesystem::temp_directory_path();
    std::filesystem::path graph_dir = tmp / "heiconnect_test_graphs";
    std::filesystem::create_directories(graph_dir);
    cycle_graph.write_to_file_metis(graph_dir / "cycle_graph.graph");
    cycle_graph.write_to_file_graphML(graph_dir / "cycle_graph.xml", map_to_original_graph);

    auto g = graph::GraphPair{};
    g.read_graph(graph_dir / "cycle_graph.graph",
                 graph_dir / "cycle_graph.xml");
    g.add_links(1, 1.f, 0); // all links = 1.0
    auto solution = solver::greedy_heuristic_strong(g);
}