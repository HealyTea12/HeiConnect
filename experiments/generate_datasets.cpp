#include <boost/program_options.hpp>

#include "HeiConnect/data_structures/graph_utils.hpp"

namespace po = boost::program_options;

void write_cycle_graphs(std::filesystem::path graph_dir, size_t start, size_t stop, size_t step)
{
    auto a = 0;
    std::vector<size_t> cycle_sizes = {};
    for (size_t n = start; n <= stop; n += step)
        cycle_sizes.push_back(n);
    for (size_t n_nodes : cycle_sizes)
    {
        auto cycle_graph = create_cycle_graph_undirected(n_nodes);
        std::unordered_map<size_t, std::vector<size_t>> map_to_original_graph{};
        for (size_t i = 0; i < n_nodes; ++i)
        {
            map_to_original_graph[i] = {i};
        }
        cycle_graph.write_to_file_graphML(graph_dir / ("cycle_" + std::to_string(n_nodes) + ".xml"), map_to_original_graph);
        cycle_graph.write_to_file_metis(graph_dir / ("cycle_" + std::to_string(n_nodes) + ".graph"));
    }
}

void write_star_graphs(std::filesystem::path graph_dir, size_t start, size_t stop, size_t step)
{
    std::vector<size_t> star_sizes = {};
    for (size_t n = start; n <= stop; n += step)
        star_sizes.push_back(n);
    for (size_t n_leaves : star_sizes)
    {
        auto star_graph = create_star_graph(n_leaves);
        std::unordered_map<size_t, std::vector<size_t>> map_to_original_graph{};
        for (size_t i = 0; i < n_leaves + 1; ++i)
        {
            map_to_original_graph[i] = {i};
        }
        star_graph.write_to_file_graphML(graph_dir / ("star_" + std::to_string(n_leaves) + ".xml"), map_to_original_graph);
        star_graph.write_to_file_metis(graph_dir / ("star_" + std::to_string(n_leaves) + ".graph"));
    }
}

int main(int argc, char **argv)
{
    po::options_description desc("Options");
    desc.add_options()("help", "help message")("output_dir,o", po::value<std::string>())("start,s", po::value<size_t>(), "start")("stop", po::value<size_t>(), "stop")("step", po::value<size_t>(), "step")("graph_type,g", po::value<std::string>()->default_value("cycle"), "graph type: cycle or star");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help") || !vm.count("output_dir") || !vm.count("start") || !vm.count("stop") || !vm.count("step"))
    {
        std::cout << desc << std::endl;
        return 1;
    }
    auto graph_type = vm["graph_type"].as<std::string>();
    if (graph_type != "cycle" && graph_type != "star")
    {
        std::cerr << "Invalid graph type: " << graph_type << ". Must be 'cycle' or 'star'." << std::endl;
        return 1;
    }
    std::cout << "Output dir: " << vm["output_dir"].as<std::string>() << "\n";
    std::filesystem::path graph_dir = vm["output_dir"].as<std::string>();
    try
    {
        std::filesystem::create_directories(graph_dir);
    }
    catch (const std::filesystem::filesystem_error &e)
    {
        std::cerr << "Error creating output directory: " << e.what() << std::endl;
        return 1;
    }
    size_t start = vm["start"].as<size_t>();
    size_t stop = vm["stop"].as<size_t>();
    size_t step = vm["step"].as<size_t>();
    if (graph_type == "cycle")
    {
        write_cycle_graphs(graph_dir, start, stop, step);
    }
    else if (graph_type == "star")
    {
        write_star_graphs(graph_dir, start, stop, step);
    }
}