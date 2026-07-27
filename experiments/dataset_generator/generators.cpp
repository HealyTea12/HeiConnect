#include "generators.hpp"

#include "HeiConnect/data_structures/graph_utils.hpp"

#include <cassert>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace dataset_generator
{
namespace
{

using ull = unsigned long long;
using WeightFunction = std::function<double(ull, ull, const WeightedCRFGraph<> &)>;

WeightFunction create_weight_function(const Config &config)
{
    if (config.distribution == DistributionType::constant)
    {
        return [weight = config.constant_weight](ull, ull, const WeightedCRFGraph<> &)
        {
            return weight;
        };
    }

    auto engine = std::make_shared<std::mt19937>(config.seed);
    if (config.distribution == DistributionType::float_uniform)
    {
        auto lower = config.float_uniform_lower;
        auto upper = config.float_uniform_upper;
        return [engine, lower, upper](ull, ull, const WeightedCRFGraph<> &) mutable
        {
            auto fraction = std::generate_canonical<double, 53>(*engine);
            return upper - fraction * (upper - lower);
        };
    }

    auto distribution = std::make_shared<std::uniform_int_distribution<long long>>(
        config.integer_uniform_lower, config.integer_uniform_upper);
    return [engine, distribution](ull, ull, const WeightedCRFGraph<> &) mutable
    {
        return static_cast<double>((*distribution)(*engine));
    };
}

void write_graph(const WeightedCRFGraph<> &graph,
                 const std::filesystem::path &output_dir,
                 const std::string &name)
{
    std::unordered_map<std::size_t, std::vector<std::size_t>> map_to_original_graph;
    for (std::size_t node = 0; node + 1 < graph.graph.vertices.size(); ++node)
    {
        map_to_original_graph[node] = {node};
    }
    graph.write_to_file_graphML(output_dir / (name + ".xml"), map_to_original_graph);
    graph.write_to_file_metis(output_dir / (name + ".graph"));
}

void generate_graph(const Config &config)
{
    std::filesystem::create_directories(config.output);
    if (config.generator == "cycle")
    {
        write_graph(create_cycle_graph_undirected(config.nodes), config.output,
                    "cycle_" + std::to_string(config.nodes));
    }
    else if (config.generator == "star")
    {
        write_graph(create_star_graph(config.nodes - 1), config.output,
                    "star_" + std::to_string(config.nodes));
    }
    else
    {
        write_graph(create_random_tree(config.nodes, config.seed), config.output,
                    "tree_" + std::to_string(config.nodes) + "_seed_" +
                        std::to_string(config.seed));
    }
}

void write_complete_links(const WeightedCRFGraph<> &graph,
                          const WeightFunction &weight_function,
                          const std::filesystem::path &file)
{
    ull node_count = graph.graph.vertices.size() - 1;
    ull total_possible_edges = (node_count * (node_count - 1)) / 2;
    // The graph stores both directions of each undirected edge.
    ull existing_edges = graph.graph.edges.size() / 2;
    ull estimated_new_edges = total_possible_edges - existing_edges;

    std::ofstream output(file);
    if (!output)
        throw std::runtime_error("could not open output file: " + file.string());
    output.rdbuf()->pubsetbuf(nullptr, 65536); // 64KB buffer
    output << node_count << ' ' << estimated_new_edges << " 1\n";

    ull edges_written = 0;
    ull progress_interval = std::max(1ULL, estimated_new_edges / 100); // Report progress every 1%
    for (ull u = 0; u < node_count; ++u)
    {
        if (u % 100 == 0)
        {
            std::cerr << "\rProcessing vertex " << u << " / " << node_count << " ("
                      << (u * 100 / node_count) << "%)" << std::flush;
        }
        for (ull v = u + 1; v < node_count; ++v) // Start from u+1 to avoid duplicates
        {
            // better not to keep assumptions about the graph being undirected when not necessary
            bool edge_in_graph = graph.is_edge(u, v) || graph.is_edge(v, u);
            if (!edge_in_graph)
            {
                double weight = weight_function(u, v, graph);
                assert(weight > 0);
                output << u + 1 << ' ' << v + 1 << ' ' << weight << '\n';
                ++edges_written;

                if (edges_written % progress_interval == 0)
                {
                    std::cerr << "\rWritten " << edges_written << " edges..." << std::flush;
                }
            }
        }
    }
    std::cerr << "\rCompleted! Written " << edges_written << " edges to " << file << "\n";
}

void generate_links(const Config &config)
{
    if (!config.output.parent_path().empty())
        std::filesystem::create_directories(config.output.parent_path());
    auto weight_function = create_weight_function(config);
    auto graph = WeightedCRFGraph<>::read_from_file(config.input_graph);
    write_complete_links(graph, weight_function, config.output);
}

} // namespace

void generate(const Config &config)
{
    if (config.output_type == OutputType::graph)
        generate_graph(config);
    else
        generate_links(config);
}

} // namespace dataset_generator
