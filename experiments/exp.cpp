#include <random>
#include <filesystem>

#include "../src/graph.hpp"
#include "set_cover/transform_single.hpp"
#include "greedy.hpp"
#include "set_cover/graph_utils.hpp"
#include "set_cover/graph.hpp"
#include "set_cover/set_cover.hpp"
#include "tools/timer.hpp"

void experiment()
{
    auto curr_dir = std::filesystem::current_path();
    auto output_dir = curr_dir / "output";
    std::filesystem::create_directories(output_dir);
    auto output = output_dir / "set_cover_greedy_v_direct_greedy_cycle_graph_6.txt";
    WeightedCRFGraph graph = create_cycle_graph(6);
    auto link_graph = generate_links(graph, [](size_t u, size_t v)
                                     { return 1.0; });

    {
        auto timer = Timer{
            "Set Cover with Greedy Parallel on Cycle Graph of size 6",
            output};
        auto sc = construct_set_cover(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            link_graph.graph.vertices,
            link_graph.graph.edges,
            link_graph.weights);
        SetCoverSolverGreedyParallel solver{sc};
        solver.solve();
        auto solution = solver.get_solution();
    }

    auto g = graph::GraphPair{};
    auto src_dir = std::filesystem::current_path().parent_path();
    g.read_graph(src_dir / "graphs" / "misc" / "cycle-6.graph",
                 src_dir / "graphs" / "misc" / "cycle-6.xml");
    g.add_links(0, 1.0, 0); // all links = 1.0
    {
        auto timer = Timer{
            "Direct Greedy on Cycle Graph of size 6",
            output};
        auto solution = solver::greedy_strong(g);
    };
}

int main()
{
    experiment();
}
