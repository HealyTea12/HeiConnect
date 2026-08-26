#include "cli.hpp"

#include <boost/program_options.hpp>

#include <iostream>
#include <stdexcept>
#include <string>

namespace dataset_generator
{
namespace
{

namespace po = boost::program_options;

OutputType parse_output_type(const std::string &value)
{
    if (value == "graph")
        return OutputType::graph;
    if (value == "links")
        return OutputType::links;
    throw po::invalid_option_value("type must be 'graph' or 'links'");
}

DistributionType parse_distribution(const std::string &value)
{
    if (value == "constant" || value == "none")
        return DistributionType::constant;
    if (value == "float_uniform" || value == "uniform")
        return DistributionType::float_uniform;
    if (value == "integer_uniform")
        return DistributionType::integer_uniform;
    throw po::invalid_option_value(
        "distribution must be 'constant', 'float_uniform', or 'integer_uniform'");
}

void validate(const Config &config, const po::variables_map &values)
{
    if (config.output_type == OutputType::graph &&
        config.generator != "cycle" && config.generator != "star" &&
        config.generator != "tree" && config.generator != "cactus")
    {
        throw po::invalid_option_value(
            "graph generator must be 'cycle', 'star', 'tree', or 'cactus'");
    }
    if (config.output_type == OutputType::links && config.generator != "complete")
    {
        throw po::invalid_option_value("links generator must be 'complete'");
    }

    if (config.output_type == OutputType::graph)
    {
        if (!values.count("nodes"))
            throw po::required_option("nodes");
        if (config.generator == "cycle" && config.nodes < 3)
            throw po::invalid_option_value("cycle graphs require at least 3 nodes");
        if (config.generator == "star" && config.nodes < 2)
            throw po::invalid_option_value("star graphs require at least 2 nodes");
        if (config.generator == "tree" && config.nodes < 2)
            throw po::invalid_option_value("tree graphs require at least 2 nodes");
        if (config.generator == "cactus")
        {
            if (!values.count("cycles"))
                throw po::required_option("cycles");
            if (!values.count("cycle-length"))
                throw po::required_option("cycle-length");
            if (config.nodes == 0)
                throw po::invalid_option_value("cactus graphs require at least 1 node");
            if (config.cycle_length < 3)
                throw po::invalid_option_value("cactus cycles require at least 3 nodes");
            if (config.cycles > 0 &&
                config.cycle_length - 1 > (config.nodes - 1) / config.cycles)
            {
                throw po::invalid_option_value(
                    "the requested cactus cycles require more than the available nodes");
            }
        }
    }
    else
    {
        if (!values.count("input_graph"))
            throw po::required_option("input_graph");
        if (!std::filesystem::is_regular_file(config.input_graph))
            throw po::invalid_option_value("input_graph must be an existing file");
        if (config.output.extension() != ".links")
            throw po::invalid_option_value("links output must have the .links extension");
    }

    if (config.distribution == DistributionType::constant &&
        config.constant_weight <= 0.0)
        throw po::invalid_option_value("constant_weight must be positive");
    if (config.distribution == DistributionType::float_uniform &&
        (config.float_uniform_lower > config.float_uniform_upper ||
         config.float_uniform_lower < 0.0 || config.float_uniform_upper <= 0.0))
    {
        throw po::invalid_option_value(
            "float_uniform bounds must contain positive values in ascending order");
    }
    if (config.distribution == DistributionType::integer_uniform &&
        (config.integer_uniform_lower <= 0 ||
         config.integer_uniform_lower > config.integer_uniform_upper))
    {
        throw po::invalid_option_value(
            "integer_uniform bounds must be positive and in ascending order");
    }
}

} // namespace

std::optional<Config> parse_config(int argc, char **argv)
{
    po::options_description options("Dataset generator options");
    // This comments are for keeping the editor from putting everything in one line
    options.add_options()
        ("help,h", "show this help message")
        ("type,t", po::value<std::string>(), "output type: graph or links")
        ("generator,g", po::value<std::string>(), "generator: cycle, star, tree, cactus, or complete")
        ("output,o", po::value<std::string>(), "graph output directory or .links output file")
        ("input_graph,i", po::value<std::string>(), "base .graph file used to generate links")
        ("nodes,n", po::value<std::size_t>(), "number of nodes in the generated graph")
        ("cycles", po::value<std::size_t>(), "number of cycles in a cactus graph")
        ("cycle-length", po::value<std::size_t>(), "length of every cactus cycle")
        ("distribution,d", po::value<std::string>()->default_value("constant"),
         "link weight distribution: constant, float_uniform, or integer_uniform")
        ("constant_weight,c", po::value<double>()->default_value(1.0), "constant link weight")
        ("float_uniform_lower", po::value<double>()->default_value(0.0), "float distribution lower bound")
        ("float_uniform_upper", po::value<double>()->default_value(1.0), "float distribution upper bound")
        ("integer_uniform_lower", po::value<long long>()->default_value(1), "integer distribution lower bound")
        ("integer_uniform_upper", po::value<long long>()->default_value(9), "integer distribution upper bound")
        ("seed", po::value<unsigned int>()->default_value(42), "random seed");

    po::variables_map values;
    po::store(po::parse_command_line(argc, argv, options), values);
    po::notify(values);

    if (values.count("help"))
    {
        std::cout << options << '\n';
        return std::nullopt;
    }
    if (!values.count("type"))
        throw po::required_option("type");
    if (!values.count("generator"))
        throw po::required_option("generator");
    if (!values.count("output"))
        throw po::required_option("output");

    Config config;
    config.output_type = parse_output_type(values["type"].as<std::string>());
    config.generator = values["generator"].as<std::string>();
    config.output = values["output"].as<std::string>();
    if (values.count("input_graph"))
        config.input_graph = values["input_graph"].as<std::string>();
    if (values.count("nodes"))
        config.nodes = values["nodes"].as<std::size_t>();
    if (values.count("cycles"))
        config.cycles = values["cycles"].as<std::size_t>();
    if (values.count("cycle-length"))
        config.cycle_length = values["cycle-length"].as<std::size_t>();
    config.distribution = parse_distribution(values["distribution"].as<std::string>());
    config.constant_weight = values["constant_weight"].as<double>();
    config.float_uniform_lower = values["float_uniform_lower"].as<double>();
    config.float_uniform_upper = values["float_uniform_upper"].as<double>();
    config.integer_uniform_lower = values["integer_uniform_lower"].as<long long>();
    config.integer_uniform_upper = values["integer_uniform_upper"].as<long long>();
    config.seed = values["seed"].as<unsigned int>();

    validate(config, values);
    return config;
}

} // namespace dataset_generator
