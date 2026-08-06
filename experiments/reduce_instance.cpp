#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/program_options.hpp>

#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/sc_reduction/transform_single_builders.hpp"
#include "HeiConnect/set_cover/file_writer.hpp"
#include "HeiConnect/set_cover/solver_ilp.hpp"

namespace po = boost::program_options;

enum class ProblemType
{
    SET_COVER,
    MINIZINC,
    ILP,
};

ProblemType parse_type(const std::string& type)
{
    if (type == "set-cover")
        return ProblemType::SET_COVER;
    if (type == "minizinc")
        return ProblemType::MINIZINC;
    if (type == "ilp")
        return ProblemType::ILP;
    throw std::invalid_argument("type must be 'set-cover', 'minizinc', or 'ilp'");
}

void validate_input_file(const std::filesystem::path& path, const std::string& extension)
{
    if (!std::filesystem::is_regular_file(path))
        throw std::invalid_argument("input file does not exist: " + path.string());
    if (path.extension() != extension)
        throw std::invalid_argument("expected a " + extension + " file: " + path.string());
}

int main(int argc, char** argv)
{
    try
    {
        po::options_description options("Connectivity augmentation reduction options");
        options.add_options()
            ("help,h", "show this help message")
            ("graph", po::value<std::string>()->required(), "input graph in GraphML .xml format")
            ("links", po::value<std::string>()->required(), "input candidate links in .links format")
            ("output,o", po::value<std::string>()->required(), "output file path")
            ("type,t", po::value<std::string>()->required(), "target problem: set-cover, minizinc, or ilp");

        po::variables_map values;
        po::store(po::parse_command_line(argc, argv, options), values);

        if (values.count("help"))
        {
            std::cout
                << "Usage: reduce_instance --graph FILE.xml --links FILE.links "
                   "--output FILE --type TYPE\n\n"
                << options << '\n';
            return 0;
        }

        po::notify(values);

        const auto graph_file = std::filesystem::path(values["graph"].as<std::string>());
        const auto links_file = std::filesystem::path(values["links"].as<std::string>());
        const auto output_file = std::filesystem::path(values["output"].as<std::string>());
        const auto type = parse_type(values["type"].as<std::string>());

        validate_input_file(graph_file, ".xml");
        validate_input_file(links_file, ".links");

        auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
        auto links = WeightedCRFGraph<>::read_from_file_links(links_file);
        if (graph.num_vertices() != links.num_vertices())
        {
            throw std::invalid_argument(
                "graph and links have different vertex counts: " + std::to_string(graph.num_vertices()) + " and "
                + std::to_string(links.num_vertices()));
        }
        auto set_cover = construct_set_cover(
            graph.graph.vertices,
            graph.graph.edges,
            graph.weights,
            links.graph.vertices,
            links.graph.edges,
            links.weights);

        if (type == ProblemType::ILP)
        {
            auto ilp = build_set_cover_ilp_model(set_cover);
            write_set_cover_ilp_model(ilp, output_file);
        }
        else
        {
            std::ofstream output(output_file);
            if (!output)
                throw std::runtime_error("could not open output file: " + output_file.string());

            const auto format = type == ProblemType::SET_COVER
                ? SetCoverWriter::Format::DEFAULT
                : SetCoverWriter::Format::MINIZINC;
            SetCoverWriter::write(set_cover, output, format);
        }
        return 0;
    }
    catch (const po::error& error)
    {
        std::cerr << "Command line error: " << error.what() << '\n';
        return 1;
    }
    catch (const GRBException& error)
    {
        std::cerr << "Gurobi error: " << error.getMessage() << '\n';
        return 1;
    }
    catch (const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
