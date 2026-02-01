#include <filesystem>
#include <iostream>

#include "boost/program_options.hpp"

#include "algorithm_registry.hpp"
#include "algorithm_runner.hpp"
#include "experiment_utils.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"

void run_experiment(const std::filesystem::path &graph_dir,
                    const std::filesystem::path &output_file,
                    Algorithms algorithm)
{
    auto runner = create_algorithm_runner(algorithm);

    for (const auto &file : std::filesystem::directory_iterator(graph_dir))
    {
        if (!file.path().filename().string().ends_with(".xml"))
            continue;

        try
        {
            log_to_file_and_stdout("Instance: " + file.path().string(), output_file);
            auto graph = WeightedCRFGraph<>::read_from_file_graphML(file.path());
            log_to_file_and_stdout("n: " + std::to_string(graph.num_vertices()), output_file);
            log_to_file_and_stdout("m: " + std::to_string(graph.num_edges()), output_file);
            log_to_file_and_stdout("Algorithm: " + algorithm_to_string(algorithm), output_file);
            log_to_file_and_stdout("d_min: " + std::to_string(graph.min_degree()), output_file);
            log_to_file_and_stdout("d_max: " + std::to_string(graph.max_degree()), output_file);
            log_to_file_and_stdout("d_avg: " + std::to_string(graph.average_degree()), output_file);
            auto memory_usage = run_isolated_and_measure_memory_usage([&]()
                                                                      { 
                runner->run(file.path());
                std::ofstream ofs{output_file.string(), std::ios::app};
                runner->print_results(std::cout);
                runner->print_results(ofs); });
            log_to_file_and_stdout("Peak memory usage (pages): " + std::to_string(memory_usage), output_file);
            log_separator(output_file);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error processing graph " << file.path() << ": " << e.what() << std::endl;
        }
    }
}

void run_experiment_file(const std::filesystem::path &graph_file,
                         Algorithms algorithm,
                         const std::filesystem::path &output_file)
{
    auto runner = create_algorithm_runner(algorithm);
    try
    {
        log_to_file_and_stdout("Instance: " + graph_file.string(), output_file);
        auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
        log_to_file_and_stdout("n: " + std::to_string(graph.num_vertices()), output_file);
        log_to_file_and_stdout("m: " + std::to_string(graph.num_edges()), output_file);
        log_to_file_and_stdout("Algorithm: " + algorithm_to_string(algorithm), output_file);
        log_to_file_and_stdout("d_min: " + std::to_string(graph.min_degree()), output_file);
        log_to_file_and_stdout("d_max: " + std::to_string(graph.max_degree()), output_file);
        log_to_file_and_stdout("d_avg: " + std::to_string(graph.average_degree()), output_file);
        runner->run(graph_file);
        runner->print_results(std::cout);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error processing graph " << graph_file << ": " << e.what() << std::endl;
    }
}

namespace po = boost::program_options;

void print_available_algorithms(std::ostream &os)
{
    os << "Available algorithms:\n";
    for (const auto &name : ALGORITHM_NAMES)
    {
        os << "  " << name << "\n";
    }
}

int main(int argc, char **argv)
{
    try
    {
        po::options_description desc("Experiment Runner Options");
        desc.add_options()                                                              //
            ("help,h", "show this help message")                                        //
            ("output_file,o", po::value<std::string>()->required(), "output file path") //
            ("algorithm,a", po::value<std::string>()->required(), "algorithm to run")   //
            ("input_dir,i", po::value<std::string>(), "input graph directory")          //
            ("input_file,f", po::value<std::string>(), "input graph file");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            print_available_algorithms(std::cout);
            return 0;
        }

        if (vm.count("input_file") && vm.count("input_dir"))
        {
            throw std::invalid_argument("Cannot specify both input_file and input_dir.");
        }
        if (!vm.count("input_file") && !vm.count("input_dir"))
        {
            throw std::invalid_argument("Must specify either input_file or input_dir.");
        }

        po::notify(vm);

        std::string graph_file{""};
        std::string graph_dir{""};

        if (vm.count("input_dir"))
            graph_dir = vm["input_dir"].as<std::string>();
        if (vm.count("input_file"))
            graph_file = vm["input_file"].as<std::string>();
        auto output_file = std::filesystem::path(vm["output_file"].as<std::string>());
        auto algorithm_str = vm["algorithm"].as<std::string>();

        auto algorithm = algorithm_from_string(algorithm_str);
        if (!graph_dir.empty())
            run_experiment(graph_dir, output_file, algorithm);
        if (!graph_file.empty())
            run_experiment_file(graph_file, algorithm, output_file);
    }
    catch (const po::error &e)
    {
        std::cerr << "Command line error: " << e.what() << "\n";
        return 1;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
