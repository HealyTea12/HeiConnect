#include "generators.hpp"

#include "HeiConnect/data_structures/graph_utils.hpp"

#include "pugixml.hpp"

#include <cassert>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace dataset_generator
{
namespace
{

using ull = unsigned long long;
using WeightFunction = std::function<double(ull, ull, const WeightedCRFGraph<> &)>;

struct GraphMetadata
{
    std::string name;
    std::string type;
    std::string value;
};

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

// Adds graph-level metadata entries to a GraphML file that has already been written.
void add_graph_metadata(const std::filesystem::path &graph_file,
                        const std::vector<GraphMetadata> &metadata)
{
    pugi::xml_document document;
    const auto parse_result = document.load_file(graph_file.c_str());
    if (!parse_result)
        throw std::runtime_error("could not read GraphML file: " + graph_file.string());

    auto graphml_node = document.child("graphml");
    auto graph_node = graphml_node.child("graph");
    if (!graphml_node || !graph_node)
        throw std::runtime_error("GraphML file does not contain a graph: " + graph_file.string());

    const auto first_graph_child = graph_node.first_child();
    for (const auto &entry : metadata)
    {
        auto key_node = graphml_node.insert_child_before("key", graph_node);
        key_node.append_attribute("id") = entry.name.c_str();
        key_node.append_attribute("for") = "graph";
        key_node.append_attribute("attr.name") = entry.name.c_str();
        key_node.append_attribute("attr.type") = entry.type.c_str();

        auto data_node = first_graph_child
                             ? graph_node.insert_child_before("data", first_graph_child)
                             : graph_node.append_child("data");
        data_node.append_attribute("key") = entry.name.c_str();
        data_node.text().set(entry.value.c_str());
    }

    if (!document.save_file(graph_file.c_str()))
        throw std::runtime_error("could not write GraphML file: " + graph_file.string());
}

// Writes the standard graph files and then enriches the GraphML file with graph-level metadata.
void write_graph_with_metadata(const WeightedCRFGraph<> &graph,
                               const std::filesystem::path &output_dir,
                               const std::string &name,
                               const std::vector<GraphMetadata> &metadata)
{
    write_graph(graph, output_dir, name);
    add_graph_metadata(output_dir / (name + ".xml"), metadata);
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
    else if (config.generator == "tree")
    {
        write_graph(create_random_tree(config.nodes, config.seed), config.output,
                    "tree_" + std::to_string(config.nodes) + "_seed_" +
                        std::to_string(config.seed));
    }
    else if (config.generator == "cactus_cycles")
    {
        const auto graph = create_random_cactus_with_cycle_count(
            config.nodes, config.cycles, config.seed);
        const size_t bridge_count = std::count(graph.weights.begin(), graph.weights.end(), 2.0) / 2;
        const std::string name =
            "cactus_cycles_n" + std::to_string(config.nodes) +
            "_q" + std::to_string(config.cycles) +
            "_seed" + std::to_string(config.seed);
        write_graph_with_metadata(graph, config.output, name, {
            {"generator", "string", "cactus_cycles"},
            {"nodes", "long", std::to_string(config.nodes)},
            {"edges", "long", std::to_string(graph.num_edges() / 2)},
            {"cycles", "long", std::to_string(config.cycles)},
            {"proper_cycles", "long", std::to_string(config.cycles - bridge_count)},
            {"bridges", "long", std::to_string(bridge_count)},
            {"seed", "long", std::to_string(config.seed)},
        });
    }
    else if (config.generator == "cactus_variable")
    {
        const auto graph = create_random_cactus_with_cycle_sizes(
            config.nodes, config.min_cycle_size, config.max_cycle_size, config.seed);
        const size_t edge_count = graph.num_edges() / 2;
        const size_t cycle_count = edge_count - (config.nodes - 1);
        const size_t bridge_count = std::count(graph.weights.begin(), graph.weights.end(), 2.0) / 2;
        const std::string name =
            "cactus_variable_n" + std::to_string(config.nodes) +
            "_min" + std::to_string(config.min_cycle_size) +
            "_max" + std::to_string(config.max_cycle_size) +
            "_seed" + std::to_string(config.seed);
        write_graph_with_metadata(graph, config.output, name, {
            {"generator", "string", "cactus_variable"},
            {"nodes", "long", std::to_string(config.nodes)},
            {"edges", "long", std::to_string(edge_count)},
            {"cycles", "long", std::to_string(cycle_count)},
            {"bridges", "long", std::to_string(bridge_count)},
            {"min_cycle_size", "long", std::to_string(config.min_cycle_size)},
            {"max_cycle_size", "long", std::to_string(config.max_cycle_size)},
            {"seed", "long", std::to_string(config.seed)},
        });
    }
    else
    {
        const size_t cycle_node_count = config.cycles * (config.cycle_length - 1);
        const size_t bridge_count = config.nodes - 1 - cycle_node_count;
        const size_t edge_count = config.nodes - 1 + config.cycles;
        const double cycle_mass = config.nodes > 1
                                      ? static_cast<double>(cycle_node_count) /
                                            static_cast<double>(config.nodes - 1)
                                      : 0.0;
        const std::string name =
            "cactus_n" + std::to_string(config.nodes) +
            "_q" + std::to_string(config.cycles) +
            "_k" + std::to_string(config.cycle_length) +
            "_seed" + std::to_string(config.seed);
        const std::vector<GraphMetadata> metadata = {
            {"generator", "string", "simple_random_cactus"},
            {"nodes", "long", std::to_string(config.nodes)},
            {"edges", "long", std::to_string(edge_count)},
            {"cycles", "long", std::to_string(config.cycles)},
            {"cycle_length", "long", std::to_string(config.cycle_length)},
            {"bridges", "long", std::to_string(bridge_count)},
            {"cycle_budget", "long", std::to_string(cycle_node_count)},
            {"cycle_mass", "double", std::to_string(cycle_mass)},
            {"seed", "long", std::to_string(config.seed)}};

        write_graph_with_metadata(
            create_random_cactus(
                config.nodes, config.cycles, config.cycle_length, config.seed),
            config.output,
            name,
            metadata);
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
