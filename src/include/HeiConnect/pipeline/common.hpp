#pragma once
#include <string>
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