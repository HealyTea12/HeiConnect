#include <boost/program_options.hpp>

#include "HeiConnect/data_structures/graph_utils.hpp"

namespace po = boost::program_options;

using ull = unsigned long long;

void drop_contained_nodes(const std::filesystem::path &graph_file_xml)
{
    auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file_xml);
    std::filesystem::path graph_file_metis = graph_file_xml.parent_path() / (graph_file_xml.stem().string() + ".graph");
    auto graph_metis = WeightedCRFGraph<>::read_from_file(graph_file_metis);
    std::unordered_map<size_t, std::vector<size_t>> mapping{};
    for (ull u = 0; u < graph_metis.graph.vertices.size() - 1; ++u)
    {
        mapping[u] = {u};
    }
    graph.write_to_file_graphML(graph_file_xml, mapping);
}

static void write_links(std::filesystem::path filename,
                        std::vector<std::tuple<ull, ull, double>> &links)
{
    ull max_node = 0;
    for (size_t i{0}; i < links.size(); ++i)
    {
        auto [u, v, w] = links[i];
        if (u > max_node)
            max_node = u;
        if (v > max_node)
            max_node = v;
    }
    std::ofstream file{filename};
    file << max_node + 1 << " " << links.size() << " " << "1" << "\n";
    for (const auto &[u, v, w] : links)
    {
        file << u + 1 << " " << v + 1 << " " << w << "\n";
    }
}

// writes only in one direction for undirected graphs
static std::vector<std::tuple<ull, ull, double>> create_links_undirected(
    const WeightedCRFGraph<> &graph,
    std::function<double(ull, ull, const WeightedCRFGraph<> &)> weight_function)
{
    std::vector<std::tuple<ull, ull, double>> links{};
    for (ull u = 0; u < graph.graph.vertices.size() - 1; ++u)
    {
        for (ull v = 0; v < graph.graph.vertices.size() - 1; ++v)
        {
            // better not to keep assumptions about the graph being undirected when not necessary, this part of
            // the code is not performance critical
            bool edge_in_graph = graph.is_edge(u, v) || graph.is_edge(v, u);
            if (u < v && !edge_in_graph)
            {
                double weight = weight_function(u, v, graph);
                assert(weight > 0);
                links.emplace_back(u, v, weight);
            }
        }
    }
    return links;
}

static void create_links_undirected_write(
    const WeightedCRFGraph<> &graph,
    std::function<double(ull, ull, const WeightedCRFGraph<> &)> weight_function,
    std::filesystem::path &file)
{
    ull n = graph.graph.vertices.size() - 1;
    ull total_possible_edges = (n * (n - 1)) / 2;
    ull existing_edges = graph.graph.edges.size();
    ull estimated_new_edges = total_possible_edges - existing_edges;

    std::cout << "Estimated new edges to generate: " << estimated_new_edges << "\n";
    std::cout << "This may take a while...\n";

    std::ofstream outfile{file};
    outfile.rdbuf()->pubsetbuf(nullptr, 65536); // 64KB buffer
    outfile << n << " " << estimated_new_edges << " " << "1" << "\n";

    ull edges_written = 0;
    ull progress_interval = std::max(1ULL, estimated_new_edges / 100); // Report progress every 1%
    if (progress_interval == 0)
        progress_interval = 1000000; // Cap at 1M edge intervals

    for (ull u = 0; u < n; ++u)
    {
        if (u % 100 == 0)
        {
            std::cerr << "\rProcessing vertex " << u << " / " << n << " (" << (u * 100 / n) << "%)" << std::flush;
        }
        for (ull v = u + 1; v < n; ++v) // Start from u+1 to avoid duplicates
        {
            // better not to keep assumptions about the graph being undirected when not necessary
            bool edge_in_graph = graph.is_edge(u, v) || graph.is_edge(v, u);
            if (u < v && !edge_in_graph)
            {
                double weight = weight_function(u, v, graph);
                assert(weight > 0);
                outfile << u + 1 << " " << v + 1 << " " << weight << "\n";
                edges_written++;

                if (edges_written % progress_interval == 0)
                {
                    std::cerr << "\rWritten " << edges_written << " edges..." << std::flush;
                }
            }
        }
    }
    outfile.flush();
    std::cerr << "\rCompleted! Written " << edges_written << " edges to " << file << "\n";
}

static void write_cycle_graphs(std::filesystem::path graph_dir, size_t start, size_t stop, size_t step,
                               std::function<double(ull, ull, const WeightedCRFGraph<> &)> weight_function)
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
        auto links = create_links_undirected(cycle_graph, weight_function);
        write_links(graph_dir / ("cycle_" + std::to_string(n_nodes) + ".links"), links);
    }
}

static void write_star_graphs(std::filesystem::path graph_dir, size_t start, size_t stop, size_t step,
                              std::function<double(ull, ull, const WeightedCRFGraph<> &)> weight_function)
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
        auto links = create_links_undirected(star_graph, weight_function);
        write_links(graph_dir / ("star_" + std::to_string(n_leaves) + ".links"), links);
    }
}

static std::array<std::string, 2> graph_types = {
    "cycle",
    "star"};
static std::array<std::string, 3> link_distributions = {
    "none",
    "constant",
    "uniform",
};

int main(int argc, char **argv)
{
    po::options_description desc("Options");
    desc.add_options()("help", "help message")                                                          // This comments are for keeping the editor from putting everything in one line
        ("output_dir,o", po::value<std::string>())                                                      //
        ("start,s", po::value<size_t>(), "start")                                                       //
        ("stop", po::value<size_t>(), "stop")                                                           //
        ("step", po::value<size_t>(), "step")                                                           //
        ("graph_type,g", po::value<std::string>()->default_value("cycle"), "graph type: cycle or star") //
        ("link_distribution,l", po::value<std::string>()->default_value("none"));

    // constant ld
    desc.add_options()("constant_weight,c", po::value<double>()->default_value(1.0), "constant link weight (only for constant link distribution)");
    // uniform ld
    desc.add_options()("uniform_min_weight", po::value<double>()->default_value(0.0), "minimum uniform link weight (only for uniform link distribution)");
    desc.add_options()("uniform_max_weight", po::value<double>()->default_value(1.0), "maximum uniform link weight (only for uniform link distribution)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help") ||
        !vm.count("output_dir"))
    {
        std::cout << desc << std::endl;
        return 1;
    }
    std::string_view ld{vm["link_distribution"].as<std::string>()};
    if (std::find(link_distributions.begin(), link_distributions.end(), ld) == link_distributions.end())
    {
        std::cerr << "Invalid link distribution: " << ld << ".\n Must be";
        for (const auto &dist : link_distributions)
        {
            std::cerr << " '" << dist << "'";
        }
        std::cerr << std::endl;
        return 1;
    }

    auto graph_type = vm["graph_type"].as<std::string>();
    if (graph_type != "cycle" && graph_type != "star" && graph_type != "fill" && graph_type != "break")
    {
        std::cerr << "Invalid graph type: " << graph_type << ".\n Must be";
        for (const auto &gt : graph_types)
        {
            std::cerr << " '" << gt << "'";
        }
        std::cerr << std::endl;
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
    std::function<double(ull, ull, const WeightedCRFGraph<> &)> weight_function;
    if (ld == "none")
    {
        // very misleading i know
        weight_function = [](ull u, ull v, const WeightedCRFGraph<> &graph)
        {
            return 1.0;
        };
    }
    else if (ld == "constant")
    {
        double constant_weight = vm["constant_weight"].as<double>();
        weight_function = [constant_weight](ull u, ull v, const WeightedCRFGraph<> &graph)
        {
            return constant_weight;
        };
    }
    else if (ld == "uniform")
    {
        double min_weight = vm["uniform_min_weight"].as<double>();
        double max_weight = vm["uniform_max_weight"].as<double>();
        if (min_weight >= max_weight)
        {
            throw po::invalid_option_value("uniform_min_weight must be less than uniform_max_weight");
        }
        if (min_weight < 0.0 || max_weight <= 0.0)
        {
            throw po::invalid_option_value("uniform link weights must be positive");
        }
        // Capture by value or use shared_ptr to avoid dangling references
        auto mt = std::make_shared<std::mt19937>(42);
        auto dist = std::make_shared<std::uniform_real_distribution<double>>(min_weight + 1e-6, max_weight);
        weight_function = [mt, dist](ull u, ull v, const WeightedCRFGraph<> &graph) mutable
        {
            return (*dist)(*mt);
        };
    }

    if (graph_type == "break")
    {
        for (auto &file : std::filesystem::directory_iterator(graph_dir))
        {
            if (file.path().extension() == ".xml")
            {
                auto graph = WeightedCRFGraph<>::read_from_file_graphML(file.path());
                graph.write_to_file_metis(file.path().parent_path() / (file.path().stem().string() + ".graph"));
                drop_contained_nodes(file.path());
            }
        }
        return 0;
    }

    else if (graph_type == "fill")
    {
        for (auto &file : std::filesystem::directory_iterator(graph_dir))
        {
            if (file.path().extension() == ".graph")
            {
                auto links_file = std::filesystem::path(file.path().parent_path() / (file.path().stem().string() + ".links"));
                std::cout << "\n=== Processing file: " << file.path() << " ===\n";
                auto graph = WeightedCRFGraph<>::read_from_file(file.path());
                ull n = graph.graph.vertices.size() - 1;
                ull m = graph.graph.edges.size();
                std::cout << "n: " << n << ", m: " << m << "\n";
                std::cout << "Density: " << (2.0 * m / (n * (n - 1))) * 100 << "%\n";
                create_links_undirected_write(graph, weight_function, links_file);
            }
        }
        return 0;
    }
    if (graph_type == "cycle")
    {
        size_t start = vm["start"].as<size_t>();
        size_t stop = vm["stop"].as<size_t>();
        size_t step = vm["step"].as<size_t>();
        write_cycle_graphs(graph_dir, start, stop, step, weight_function);
    }
    else if (graph_type == "star")
    {
        size_t start = vm["start"].as<size_t>();
        size_t stop = vm["stop"].as<size_t>();
        size_t step = vm["step"].as<size_t>();
        write_star_graphs(graph_dir, start, stop, step, weight_function);
    }
}