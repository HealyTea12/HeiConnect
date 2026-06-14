#pragma once

#include <vector>
#include <any>
#include <optional>
#include <string>
#include <functional>
#include <chrono>

#include "HeiConnect/pipeline/common.hpp"

class RuntimePipeline
{
public:
    using TimePoint = std::chrono::time_point<std::chrono::steady_clock>;
    using Duration = std::chrono::duration<double>;
    using Clock = std::chrono::steady_clock;

    virtual ~RuntimePipeline() = default;
    RuntimePipeline& add_step(std::any step)
    {
        m_pipeline.push_back(step);
        return *this;
    }
    std::any run(std::any input)
    {
        std::any current = input;
        for (const auto& step : m_pipeline)
        {
            TimePoint start = Clock::now();
            current = std::any_cast<std::function<std::any(std::any)>>(step.second)(current);
            TimePoint end = Clock::now();
            Duration duration = end - start;
        }
    }

private:
    std::vector<std::pair<std::optional<std::string>, std::function<std::any(std::any)>>> m_pipeline;
};


#pragma once

#include "pipeline.hpp"

#include <any>
#include <chrono>
#include <functional>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

namespace runtime_pipeline_detail
{
    template<class>
    inline constexpr bool always_false_v = false;

    template<class T>
    struct function_traits
    {
        static constexpr bool is_supported = false;
    };

    template<class R, class... Args>
    struct function_traits<R(Args...)>
    {
        static constexpr bool is_supported = true;
        using return_type = R;
        using args_tuple = std::tuple<Args...>;
    };

    template<class R, class... Args>
    struct function_traits<R (*)(Args...)> : function_traits<R(Args...)>
    {};

    template<class R, class... Args>
    struct function_traits<R (*)(Args...) noexcept> : function_traits<R(Args...)>
    {};

    template<class C, class R, class... Args>
    struct function_traits<R (C::*)(Args...)> : function_traits<R(Args...)>
    {};

    template<class C, class R, class... Args>
    struct function_traits<R (C::*)(Args...) const> : function_traits<R(Args...)>
    {};

    template<class C, class R, class... Args>
    struct function_traits<R (C::*)(Args...) noexcept> : function_traits<R(Args...)>
    {};

    template<class C, class R, class... Args>
    struct function_traits<R (C::*)(Args...) const noexcept> : function_traits<R(Args...)>
    {};

    template<class R, class... Args>
    struct function_traits<std::function<R(Args...)>> : function_traits<R(Args...)>
    {};

    template<class F>
        requires requires { &std::remove_cvref_t<F>::operator(); }
    struct function_traits<F> : function_traits<decltype(&std::remove_cvref_t<F>::operator())>
    {};

    template<class T>
    std::string type_name()
    {
        return typeid(std::remove_cvref_t<T>).name();
    }

    inline std::string type_name(const std::any& value)
    {
        return value.has_value() ? value.type().name() : std::string{"<empty std::any>"};
    }

    template<class Expected>
    [[noreturn]] void throw_bad_state_type(const std::any& state)
    {
        std::ostringstream oss;
        oss << "RuntimePipeline stage expected state of type " << type_name<Expected>() << " but received "
            << type_name(state);
        throw std::invalid_argument(oss.str());
    }

    template<class... Args>
    using decayed_tuple_t = std::tuple<std::remove_cvref_t<Args>...>;

    template<class Arg>
    decltype(auto) any_cast_arg(std::any& state)
    {
        using Value = std::remove_cvref_t<Arg>;

        if constexpr (std::is_same_v<Value, std::any>)
        {
            if constexpr (std::is_lvalue_reference_v<Arg>)
            {
                return static_cast<Arg>(state);
            }
            else
            {
                return std::move(state);
            }
        }
        else
        {
            auto* value = std::any_cast<Value>(&state);
            if (value == nullptr)
            {
                throw_bad_state_type<Value>(state);
            }

            if constexpr (std::is_lvalue_reference_v<Arg>)
            {
                return static_cast<Arg>(*value);
            }
            else
            {
                return static_cast<Value&&>(*value);
            }
        }
    }

    template<class Stage, class... Args>
    decltype(auto) invoke_from_tuple(Stage& stage, std::any& state)
    {
        using Tuple = decayed_tuple_t<Args...>;
        auto* tuple = std::any_cast<Tuple>(&state);
        if (tuple == nullptr)
        {
            throw_bad_state_type<Tuple>(state);
        }

        return std::apply(
            [&stage](auto&&... values) -> decltype(auto) {
                return std::invoke(stage, std::forward<decltype(values)>(values)...);
            },
            std::move(*tuple));
    }

    template<class... Args>
    struct first_arg;

    template<class First, class... Rest>
    struct first_arg<First, Rest...>
    {
        using type = First;
    };

    template<class Stage, class... Args>
    decltype(auto) invoke_stage(Stage& stage, std::any& state)
    {
        if constexpr (sizeof...(Args) == 0)
        {
            using Tuple = std::tuple<>;
            if (std::any_cast<Tuple>(&state) == nullptr)
            {
                throw_bad_state_type<Tuple>(state);
            }
            return std::invoke(stage);
        }
        else if constexpr (sizeof...(Args) == 1)
        {
            using Tuple = decayed_tuple_t<Args...>;
            if (std::any_cast<Tuple>(&state) != nullptr)
            {
                return invoke_from_tuple<Stage, Args...>(stage, state);
            }

            using Arg = typename first_arg<Args...>::type;
            return std::invoke(stage, any_cast_arg<Arg>(state));
        }
        else
        {
            return invoke_from_tuple<Stage, Args...>(stage, state);
        }
    }

    template<class Emitted>
    std::optional<StageMetrics> normalize_metrics(Emitted&& emitted)
    {
        using MetricsType = std::remove_cvref_t<Emitted>;

        if constexpr (std::is_same_v<MetricsType, std::optional<StageMetrics>>)
        {
            if (emitted.has_value())
            {
                return std::move(*emitted);
            }
            return std::nullopt;
        }
        else if constexpr (std::is_same_v<MetricsType, StageMetrics>)
        {
            return std::forward<Emitted>(emitted);
        }
        else
        {
            static_assert(
                always_false_v<MetricsType>,
                "emit_metrics() must return StageMetrics or std::optional<StageMetrics>.");
        }
    }

    template<class Stage, class Result>
    std::optional<StageMetrics> emit_metrics_if_available(Stage& stage, Result& result)
    {
        if constexpr (requires { stage.emit_metrics(); })
        {
            auto emitted_metrics = stage.emit_metrics();
            return normalize_metrics(std::move(emitted_metrics));
        }
        else if constexpr (requires { stage.emit_metrics(result); })
        {
            auto emitted_metrics = stage.emit_metrics(result);
            return normalize_metrics(std::move(emitted_metrics));
        }
        else
        {
            return std::nullopt;
        }
    }

    template<class Stage>
    std::string stage_name(const Stage& stage)
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

    inline std::any into_any(std::any value)
    {
        return value;
    }

    template<class T>
    std::any into_any(T&& value)
    {
        return std::any{std::forward<T>(value)};
    }
}

class RuntimePipeline
{
public:
    using TimePoint = std::chrono::time_point<std::chrono::steady_clock>;
    using Duration = std::chrono::duration<double>;
    using Clock = std::chrono::steady_clock;

    using PerStageMetrics = ::PerStageMetrics;
    using PipelineMetrics = ::PipelineMetrics;

    virtual ~RuntimePipeline() = default;

    RuntimePipeline& add_step(std::function<std::any(std::any)> step)
    {
        return add_raw_step(std::nullopt, std::move(step));
    }

    RuntimePipeline& add_step(std::string name, std::function<std::any(std::any)> step)
    {
        return add_raw_step(std::move(name), std::move(step));
    }

    RuntimePipeline& add_step(const char* name, std::function<std::any(std::any)> step)
    {
        return add_raw_step(std::string{name}, std::move(step));
    }

    RuntimePipeline& add_step(std::any step)
    {
        if (auto* function = std::any_cast<std::function<std::any(std::any)>>(&step))
        {
            return add_raw_step(std::nullopt, *function);
        }

        using NamedStep = std::pair<std::optional<std::string>, std::function<std::any(std::any)>>;
        if (auto* named_step = std::any_cast<NamedStep>(&step))
        {
            return add_raw_step(named_step->first, named_step->second);
        }

        throw std::invalid_argument(
            "RuntimePipeline::add_step(std::any) expects either "
            "std::function<std::any(std::any)> or "
            "std::pair<std::optional<std::string>, std::function<std::any(std::any)>>.");
    }

    template<class Stage>
    RuntimePipeline& add_step(Stage stage)
    {
        return add_typed_step(std::nullopt, std::move(stage));
    }

    template<class Stage>
    RuntimePipeline& add_step(std::string name, Stage stage)
    {
        return add_typed_step(std::move(name), std::move(stage));
    }

    template<class Stage>
    RuntimePipeline& add_step(const char* name, Stage stage)
    {
        return add_typed_step(std::string{name}, std::move(stage));
    }

    template<class... Args, class Stage>
    RuntimePipeline& add_step_as(Stage stage)
    {
        return add_typed_step_with_signature<Stage, Args...>(std::nullopt, std::move(stage));
    }

    template<class... Args, class Stage>
    RuntimePipeline& add_step_as(std::string name, Stage stage)
    {
        return add_typed_step_with_signature<Stage, Args...>(std::move(name), std::move(stage));
    }

    template<class... Args, class Stage>
    RuntimePipeline& add_step_as(const char* name, Stage stage)
    {
        return add_typed_step_with_signature<Stage, Args...>(std::string{name}, std::move(stage));
    }

    std::pair<std::any, PipelineMetrics> run(std::any input)
    {
        return run_any(std::move(input));
    }

    template<class... InitialArgs>
    std::pair<std::any, PipelineMetrics> run(InitialArgs&&... args)
    {
        auto initial_state = std::make_tuple(std::forward<InitialArgs>(args)...);
        return run_any(std::any{std::move(initial_state)});
    }

    static void print_pipeline_metrics(const PipelineMetrics& metrics, std::ostream& os)
    {
        os << "Pipeline metrics:\n";
        for (size_t i = 0; i < metrics.stages.size(); ++i)
        {
            const auto& m = metrics.stages[i];
            os << " Stage " << i << " (" << m.stage_name << "): " << m.duration_seconds << "s";
            if (m.metrics.has_value())
            {
                for (const auto& metric : *m.metrics)
                {
                    os << ", " << metric.name << "=" << metric.printable_value;
                }
            }
            os << "\n";
        }
    }

private:
    struct StageResult
    {
        std::any state;
        std::optional<StageMetrics> metrics;
    };

    struct StageEntry
    {
        std::string name;
        std::function<StageResult(std::any)> run;
    };

    std::vector<StageEntry> m_pipeline;

    RuntimePipeline& add_raw_step(std::optional<std::string> name, std::function<std::any(std::any)> step)
    {
        StageEntry entry;
        entry.name = name.value_or(std::string{typeid(step).name()});
        entry.run = [step = std::move(step)](std::any current) mutable -> StageResult {
            return StageResult{step(std::move(current)), std::nullopt};
        };

        m_pipeline.emplace_back(std::move(entry));
        return *this;
    }

    template<class Stage>
    RuntimePipeline& add_typed_step(std::optional<std::string> name, Stage stage)
    {
        using StageType = std::remove_cvref_t<Stage>;
        using Traits = runtime_pipeline_detail::function_traits<StageType>;

        if constexpr (!Traits::is_supported)
        {
            static_assert(
                runtime_pipeline_detail::always_false_v<StageType>,
                "RuntimePipeline::add_step() could not infer this callable's argument types. "
                "Use add_step_as<Arg1, Arg2, ...>(stage) for generic or overloaded callables.");
        }
        else
        {
            return add_typed_step_from_tuple<StageType>(
                std::move(name),
                std::move(stage),
                typename Traits::args_tuple{});
        }
    }

    template<class Stage, class... Args>
    RuntimePipeline& add_typed_step_from_tuple(std::optional<std::string> name, Stage stage, std::tuple<Args...>)
    {
        return add_typed_step_with_signature<Stage, Args...>(std::move(name), std::move(stage));
    }

    template<class Stage, class... Args>
    RuntimePipeline& add_typed_step_with_signature(std::optional<std::string> name, Stage stage)
    {
        using StageType = std::remove_cvref_t<Stage>;
        using InvokeResult = std::invoke_result_t<StageType&, Args...>;

        if constexpr (std::is_void_v<InvokeResult>)
        {
            static_assert(
                runtime_pipeline_detail::always_false_v<StageType>,
                "RuntimePipeline stages must return the next pipeline state, just like Pipeline stages.");
        }
        else
        {
            auto stage_ptr = std::make_shared<StageType>(std::move(stage));

            StageEntry entry;
            entry.name = name.value_or(runtime_pipeline_detail::stage_name(*stage_ptr));
            entry.run = [stage_ptr](std::any current) mutable -> StageResult {
                auto result = runtime_pipeline_detail::invoke_stage<StageType, Args...>(*stage_ptr, current);
                auto stage_metrics = runtime_pipeline_detail::emit_metrics_if_available(*stage_ptr, result);

                return StageResult{runtime_pipeline_detail::into_any(std::move(result)), std::move(stage_metrics)};
            };

            m_pipeline.emplace_back(std::move(entry));
            return *this;
        }
    }

    std::pair<std::any, PipelineMetrics> run_any(std::any initial_state)
    {
        PipelineMetrics metrics;
        std::any current = std::move(initial_state);

        for (auto& stage : m_pipeline)
        {
            TimePoint start = Clock::now();
            StageResult result = stage.run(std::move(current));
            TimePoint end = Clock::now();

            PerStageMetrics psm;
            psm.duration_seconds = Duration{end - start}.count();
            psm.metrics = std::move(result.metrics);
            psm.stage_name = stage.name;
            metrics.stages.emplace_back(std::move(psm));

            current = std::move(result.state);
        }

        return std::pair{std::move(current), std::move(metrics)};
    }
};
