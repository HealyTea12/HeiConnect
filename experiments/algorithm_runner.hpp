#pragma once

#include <filesystem>
#include <memory>
#include "algorithm_registry.hpp"
#include "HeiConnect/pipeline/pipeline.hpp"

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
        std::optional<double> time_data_reduction;
    } result;
    Pipeline<>::PipelineMetrics pipeline_metrics;

    virtual ~AlgorithmRunner() = default;
    virtual void run(const std::filesystem::path& graph_file) = 0;
    virtual void print_results(std::ostream& os)
    {
        // If pipeline metrics exist, use the pipeline printer. Otherwise, do nothing.
        if (!pipeline_metrics.stages.empty())
        {
            Pipeline<>::print_pipeline_metrics(pipeline_metrics, os);
        }
    }
};

std::unique_ptr<AlgorithmRunner> create_algorithm_runner(Algorithms algorithm);
