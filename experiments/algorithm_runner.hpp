#pragma once

#include <filesystem>
#include <memory>
#include <optional>
#include <ostream>
#include <tuple>
#include <vector>

#include "HeiConnect/pipeline/pipeline.hpp"

class AlgorithmRunner
{
public:
    using Solution = std::tuple<std::vector<size_t>, double>;
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
        std::optional<double> time_data_reduction;
    } result;
    Pipeline<>::PipelineMetrics pipeline_metrics;

    virtual ~AlgorithmRunner() = default;
    virtual void run(
        const std::filesystem::path& graph_file,
        const std::filesystem::path& link_file,
        const std::filesystem::path& output_dir) = 0;
    virtual void print_results(std::ostream& os)
    {
        // If pipeline metrics exist, use the pipeline printer. Otherwise, print the results in a simple format.
        if (!pipeline_metrics.stages.empty())
        {
            Pipeline<>::print_pipeline_metrics(pipeline_metrics, os);
        }
        else
        {
            // Print results in a simple format
            if (result.solution_cost.has_value())
            {
                os << "solution.cost=" << *result.solution_cost << "\n";
            }
            if (result.solution_size.has_value())
            {
                os << "solution.size=" << *result.solution_size << "\n";
                if (result.time_total.has_value())
                {
                    os << "algorithm.total_time_seconds=" << *result.time_total << "\n";
                }
            }
        }
    }
};
