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
    const std::filesystem::path& result_file_name,
    bool log_stdout,
    const ParamMap& params)
{
    auto runner = global_registry.create(algorithm, params);
    const auto result_file = output_dir / result_file_name;
    try
    {
        auto memory_usage = run_isolated_and_measure_memory_usage([&]() {
            const auto start = std::chrono::steady_clock::now();
            runner->run(graph_file, link_file, output_dir);
            const auto total_time = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            std::ofstream ofs{result_file.string(), std::ios::app};
            if (log_stdout)
            {
                std::cout << "run.instance=" << graph_file.string() << "\n\n";
                runner->print_results(std::cout);
                std::cout << "run.total_time_seconds=" << total_time << "\n";
            }
            ofs << "run.instance=" << graph_file.string() << "\n\n";
            runner->print_results(ofs);
            ofs << "run.total_time_seconds=" << total_time << "\n";
        });
        log_to_file_and_stdout("run.peak_memory_pages=" + std::to_string(memory_usage), result_file);
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
            ("result_file_name",
             po::value<std::string>()->default_value("res.txt"),
             "name of the result file inside the output directory") //
            ("algorithm,a", po::value<std::string>()->required(), "algorithm to run") //
            ("input_file,f",
             po::value<std::vector<std::string>>()->multitoken()->required(),
             "input graph file and links file")(
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

        po::notify(vm);

        const auto input_files = vm["input_file"].as<std::vector<std::string>>();
        if (input_files.size() != 2)
        {
            throw std::invalid_argument("input_file requires a graph file and a links file.");
        }
        ParamMap params;

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
        const auto result_file_name = std::filesystem::path(vm["result_file_name"].as<std::string>());
        if (result_file_name.empty() || result_file_name.has_parent_path())
        {
            throw std::invalid_argument("result_file_name must be a file name without a directory");
        }
        std::filesystem::create_directories(output_dir);
        auto algorithm_str = vm["algorithm"].as<std::string>();
        const auto names = global_registry.names();
        if (std::find(names.begin(), names.end(), algorithm_str) == names.end())
        {
            std::cerr << "Error: Algorithm '" << algorithm_str << "' not found in registry.\n";
            print_available_algorithms(std::cerr);
            return 1;
        }

        run_experiment_file(
            input_files[0],
            input_files[1],
            algorithm_str,
            output_dir,
            result_file_name,
            vm["log_stdout"].as<bool>(),
            params);
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
