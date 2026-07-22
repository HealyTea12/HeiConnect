#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <vector>

#include "boost/program_options.hpp"

#include "algorithm_registry.hpp"
#include "algorithm_runner.hpp"
#include "experiment_utils.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"

void run_experiment_file(
    const std::filesystem::path& graph_file,
    const std::filesystem::path& link_file,
    const std::string_view algorithm,
    const std::filesystem::path& output_dir,
    bool log_stdout,
    const ParamMap& params);

void run_experiment(
    const std::filesystem::path& graph_dir,
    const std::filesystem::path& output_dir,
    const std::string_view algorithm,
    bool log_stdout,
    const ParamMap& params)
{
    for (const auto& file : std::filesystem::directory_iterator(graph_dir))
    {
        if (!file.path().filename().string().ends_with(".xml"))
            continue;

        try
        {
            const auto link_file = file.path().parent_path() / (file.path().stem().string() + ".links");
            run_experiment_file(file.path(), link_file, algorithm, output_dir, log_stdout, params);
            log_separator(output_dir / "res.txt");
        }
        catch (const std::exception& e)
        {
            log_separator(output_dir / "res.txt");
            std::cerr << "Error processing graph " << file.path() << ": " << e.what() << std::endl;
        }
    }
}

void run_experiment_file(
    const std::filesystem::path& graph_file,
    const std::filesystem::path& link_file,
    const std::string_view algorithm,
    const std::filesystem::path& output_dir,
    bool log_stdout,
    const ParamMap& params)
{
    auto runner = global_registry.create(algorithm, params);
    const auto result_file = output_dir / "res.txt";
    try
    {
        auto memory_usage = run_isolated_and_measure_memory_usage([&]() {
            const auto start = std::chrono::steady_clock::now();
            runner->run(graph_file, link_file, output_dir);
            const auto total_time = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            std::ofstream ofs{result_file.string(), std::ios::app};
            if (log_stdout)
            {
                std::cout << "Instance: " << graph_file.string() << "\n";
                runner->print_results(std::cout);
                std::cout << "Total time: " << total_time << "s\n";
            }
            ofs << "Instance: " << graph_file.string() << "\n";
            runner->print_results(ofs);
            ofs << "Total time: " << total_time << "s\n";
        });
        log_to_file_and_stdout("Peak memory usage (pages): " + std::to_string(memory_usage), result_file);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error processing graph " << graph_file << ": " << e.what() << std::endl;
    }
}

namespace po = boost::program_options;

void print_available_algorithms(std::ostream& os)
{
    os << "Available algorithms:\n";
    for (const auto& name : global_registry.names())
    {
        os << "  " << name << "\n";
    }
}

int main(int argc, char** argv)
{
    try
    {
        po::options_description desc("Experiment Runner Options");
        desc.add_options() //
            ("help,h", "show this help message") //
            ("output_dir,o", po::value<std::string>()->required(), "output directory path") //
            ("algorithm,a", po::value<std::string>()->required(), "algorithm to run") //
            ("input_dir,i", po::value<std::string>(), "input graph directory") //
            ("input_file,f", po::value<std::vector<std::string>>()->multitoken(), "input graph file and links file")(
                "param,p",
                po::value<std::vector<std::string>>()->multitoken(),
                "algorithm parameters as key=value")(
                "log_stdout",
                po::bool_switch()->default_value(false),
                "whether to also log results to stdout");

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

        std::vector<std::string> input_files;
        std::string graph_dir{""};
        ParamMap params;

        if (vm.count("input_dir"))
            graph_dir = vm["input_dir"].as<std::string>();
        if (vm.count("input_file"))
        {
            input_files = vm["input_file"].as<std::vector<std::string>>();
            if (input_files.size() != 2)
                throw std::invalid_argument("input_file requires a graph file and a links file.");
        }
        if (vm.count("param"))
        {
            for (const auto& parameter : vm["param"].as<std::vector<std::string>>())
            {
                const auto separator = parameter.find('=');
                if (separator == std::string::npos)
                    throw std::invalid_argument("Algorithm parameter must have the form key=value: " + parameter);

                params[parameter.substr(0, separator)] = parameter.substr(separator + 1);
            }
        }
        const auto output_dir = std::filesystem::path(vm["output_dir"].as<std::string>());
        std::filesystem::create_directories(output_dir);
        auto algorithm_str = vm["algorithm"].as<std::string>();
        const auto names = global_registry.names();
        if (std::find(names.begin(), names.end(), algorithm_str) == names.end())
        {
            std::cerr << "Error: Algorithm '" << algorithm_str << "' not found in registry.\n";
            print_available_algorithms(std::cerr);
            return 1;
        }

        if (!graph_dir.empty())
            run_experiment(graph_dir, output_dir, algorithm_str, vm["log_stdout"].as<bool>(), params);
        if (!input_files.empty())
            run_experiment_file(
                input_files[0], input_files[1], algorithm_str, output_dir, vm["log_stdout"].as<bool>(), params);
    }
    catch (const po::error& e)
    {
        std::cerr << "Command line error: " << e.what() << "\n";
        return 1;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
