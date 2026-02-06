#pragma once

#include <filesystem>
#include <memory>
#include <optional>
#include "algorithm_registry.hpp"

class AlgorithmRunner
{
public:
    struct Result
    {
        std::optional<double> solution_cost;
        std::optional<double> solution_cost_trimmed;
        std::optional<double> solution_cost_ls;
        std::optional<size_t> solution_size;
        std::optional<size_t> solution_size_trimmed;
        std::optional<size_t> solution_size_ls;
        std::optional<double> time_reduction;
        std::optional<double> time_solving;
        std::optional<double> time_trimming;
        std::optional<double> time_ls;
        std::optional<double> time_total;
    } result;

    virtual ~AlgorithmRunner() = default;
    virtual void run(const std::filesystem::path &graph_file) = 0;
    virtual void print_results(std::ostream &os)
    {
        auto print_opt = [&os](const std::string &label, const auto &value)
        {
            if (value.has_value())
                os << label << ": " << *value << "\n";
        };
        print_opt("Solution cost", result.solution_cost);
        print_opt("Solution cost (trimmed)", result.solution_cost_trimmed);
        print_opt("Solution cost (local search)", result.solution_cost_ls);
        print_opt("Solution size", result.solution_size);
        print_opt("Solution size (trimmed)", result.solution_size_trimmed);
        print_opt("Solution size (local search)", result.solution_size_ls);
        print_opt("Time reduction (s)", result.time_reduction);
        print_opt("Time solving (s)", result.time_solving);
        print_opt("Time total (s)", result.time_total);
        print_opt("Time trimming (s)", result.time_trimming);
        print_opt("Time local search (s)", result.time_ls);
    }
};

std::unique_ptr<AlgorithmRunner> create_algorithm_runner(Algorithms algorithm);
