#pragma once

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <tuple>

#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/pipeline/common.hpp"
#include "HeiConnect/set_cover/common.hpp"
#include "HeiConnect/set_cover/set_cover_cyc.hpp"
#include "HeiConnect/set_cover/set_cover_double.hpp"
#include "HeiConnect/set_cover/solver_greedy_context.hpp"

namespace
{
    template<typename SetCoverType, typename ContextType>
    auto make_bound_context_stage()
    {
        return [](std::shared_ptr<const SetCoverType> set_cover) {
            ContextType context{*set_cover};
            BoundContext bound_context{std::move(context), USSolution{}};
            return std::tuple{std::move(set_cover), std::move(bound_context)};
        };
    }

    template<typename Solver>
    auto make_solver_stage(Solver& solver)
    {
        return [&solver](auto set_cover, auto context) {
            solver.solve(*set_cover, context);
            return std::tuple{std::move(set_cover), std::move(context)};
        };
    }
} // namespace

class GraphMetricsCalculator
{
public:
    static constexpr std::string_view name = "Graph Metrics";

    auto run(const WeightedCRFGraph<>& graph, const WeightedCRFGraph<>& link_graph)
    {
        m_metrics = StageMetrics{
            {"n", std::to_string(graph.num_vertices())},
            {"m", std::to_string(graph.num_edges())},
            {"d_min", std::to_string(graph.min_degree())},
            {"d_max", std::to_string(graph.max_degree())},
            {"d_avg", std::to_string(graph.average_degree())},
            {"num_links", std::to_string(link_graph.num_edges())},
        };
        return std::tuple{graph, link_graph};
    }
    auto operator()(const WeightedCRFGraph<>& graph, const WeightedCRFGraph<>& link_graph)
    {
        return run(graph, link_graph);
    }

    std::optional<StageMetrics> emit_metrics() const
    {
        if (!m_metrics.has_value())
        {
            return std::nullopt;
        }
        else
        {
            return m_metrics;
        }
    }

private:
    std::optional<StageMetrics> m_metrics;
};

class ReadFromFileCSRStage
{
public:
    static constexpr std::string_view name = "Read from file (CSR)";

    std::tuple<WeightedCRFGraph<>, WeightedCRFGraph<>> operator()(
        const std::filesystem::path& graph_file,
        const std::filesystem::path& link_file)
    {
        auto graph = WeightedCRFGraph<>::read_from_file_graphML(graph_file);
        auto link_graph = WeightedCRFGraph<>::read_from_file_links(link_file);
        return std::make_tuple(std::move(graph), std::move(link_graph));
    }
};
