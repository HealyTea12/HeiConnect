#pragma once
#include <functional>
#include <tuple>
#include <utility>
#include <optional>
#include <string>
#include <string_view>
#include <iostream>
#include <vector>
#include <type_traits>
#include <typeinfo>
#include <omp.h>

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

template<class... Stages>
class Pipeline
{
public:
    explicit Pipeline(Stages... stages) : m_stages(std::move(stages)...)
    {}

    using PerStageMetrics = ::PerStageMetrics;
    using PipelineMetrics = ::PipelineMetrics;

    template<class... InitialArgs>
    auto run(InitialArgs&&... args)
    {
        PipelineMetrics metrics;
        auto initial_state = std::make_tuple(std::forward<InitialArgs>(args)...);
        auto final_state = run_impl<0>(std::move(initial_state), metrics);
        return std::pair{std::move(final_state), std::move(metrics)};
    }

    static void print_pipeline_metrics(const PipelineMetrics &metrics, std::ostream &os)
    {
        os << "Pipeline metrics:\n";
        for (size_t i = 0; i < metrics.stages.size(); ++i)
        {
            const auto &m = metrics.stages[i];
            os << " Stage " << i << " (" << m.stage_name << "): " << m.duration_seconds << "s";
            if (m.metrics.has_value())
            {
                for (const auto &metric : *m.metrics)
                {
                    os << ", " << metric.name << "=" << metric.printable_value;
                }
            }
            os << "\n";
        }
    }

private:
    std::tuple<Stages...> m_stages;

    template<typename Stage, typename State>
    static auto invoke_stage(Stage &stage, State &&state)
    {
        if constexpr (requires { typename std::tuple_size<std::remove_cvref_t<State>>::type; })
        {
            return std::apply(stage, std::forward<State>(state));
        }
        else
        {
            return std::invoke(stage, std::forward<State>(state));
        }
    }

    template<typename Stage>
    static auto stage_name(const Stage &stage)
    {
        using StageType = std::remove_cvref_t<Stage>;
        if constexpr (requires { StageType::name; })
        {
            return std::string{StageType::name};
        }
        else
        {
            return std::string{typeid(stage).name()};
        }
    }

    template<std::size_t I, class State>
    auto run_impl(State state, PipelineMetrics &metrics)
    {
        if constexpr (I == sizeof...(Stages))
        {
            return state;
        }
        else
        {
            auto &stage = std::get<I>(m_stages);
            double start = omp_get_wtime();
            auto next_state = invoke_stage(stage, std::move(state));
            double end = omp_get_wtime();

            std::optional<StageMetrics> stage_metrics;
            if constexpr (requires { stage.emit_metrics(); })
            {
                auto emitted_metrics = stage.emit_metrics();
                if (emitted_metrics.has_value())
                {
                    stage_metrics = std::move(*emitted_metrics);
                }
            }
            else if constexpr (requires { stage.emit_metrics(next_state); })
            {
                auto emitted_metrics = stage.emit_metrics(next_state);
                if (emitted_metrics.has_value())
                {
                    stage_metrics = std::move(*emitted_metrics);
                }
            }

            PerStageMetrics psm;
            psm.duration_seconds = end - start;
            psm.metrics = std::move(stage_metrics);
            psm.stage_name = stage_name(stage);
            metrics.stages.emplace_back(std::move(psm));

            return run_impl<I + 1>(std::move(next_state), metrics);
        }
    }
};

template<class... Stages>
Pipeline(Stages...) -> Pipeline<Stages...>;