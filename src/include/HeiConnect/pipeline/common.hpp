#pragma once
#include <cctype>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <optional>

struct StageMetric
{
    std::string name;
    std::string printable_value;
};

using StageMetrics = std::vector<StageMetric>;

struct PerStageMetrics
{
    double duration_seconds{0.0};
    std::optional<StageMetrics> metrics;
    std::string stage_name{};
};

struct PipelineMetrics
{
    std::vector<PerStageMetrics> stages;
};

inline std::string metric_key_component(const std::string& value)
{
    std::string result;
    bool previous_was_separator = false;
    for (const unsigned char c : value)
    {
        if (std::isalnum(c))
        {
            result.push_back(static_cast<char>(std::tolower(c)));
            previous_was_separator = false;
        }
        else if (!result.empty() && !previous_was_separator)
        {
            result.push_back('_');
            previous_was_separator = true;
        }
    }
    if (!result.empty() && result.back() == '_')
    {
        result.pop_back();
    }
    return result;
}

inline void print_pipeline_metrics(const PipelineMetrics& metrics, std::ostream& os)
{
    std::unordered_map<std::string, size_t> prefix_counts;
    for (const auto& stage : metrics.stages)
    {
        const std::string base_prefix = metric_key_component(stage.stage_name);
        const size_t prefix_count = ++prefix_counts[base_prefix];
        const std::string prefix =
            prefix_count == 1 ? base_prefix : base_prefix + "_" + std::to_string(prefix_count);
        os << prefix << ".duration_seconds=" << stage.duration_seconds << "\n";
        if (stage.metrics.has_value())
        {
            for (const auto& metric : *stage.metrics)
            {
                std::string name = metric_key_component(metric.name);
                std::string value = metric.printable_value;
                if (value.size() > 1 && value.back() == 's')
                {
                    value.pop_back();
                    if (!name.ends_with("_seconds"))
                    {
                        name += "_seconds";
                    }
                }
                os << prefix << "." << name << "=" << value << "\n";
            }
        }
        os << "\n";
    }
}
