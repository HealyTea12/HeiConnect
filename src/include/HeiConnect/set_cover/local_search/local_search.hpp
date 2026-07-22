#pragma once

#include <chrono>
#include <iostream>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <variant>
#include <vector>

#include "HeiConnect/pipeline/common.hpp"
#include "HeiConnect/set_cover/solver_greedy_context.hpp"
#include "HeiConnect/set_cover/util.hpp"
#include "HeiConnect/pipeline/pipeline.hpp"


// TODO: define concept with requirements so it is not duck typed


class DefaultEvaluator
{
public:
    DefaultEvaluator(double epsilon = 1e-6) : m_epsilon(epsilon)
    {}

    template<typename SetCoverType, typename SolutionType>
    bool evaluate(const SetCoverType& set_cover, const SolutionType& solution, const SolutionType& candidate)
    {
        using SetCost = typename SetCoverType::SetCost;
        SetCost candidate_cost = HeiConnect::sc::cost(set_cover, candidate);
        SetCost solution_cost = HeiConnect::sc::cost(set_cover, solution);
        return solution_cost - candidate_cost > m_epsilon;
    }

private:
    double m_epsilon = 1e-6;
};

template<typename DestroyerType, typename RepairerType, size_t RecordMetricsLevel = 0, typename EvaluatorType = DefaultEvaluator, bool Debug = false>
class BreakAndRepairSearch
{
public:
    static constexpr std::string_view name = "Local search";

    BreakAndRepairSearch(
        DestroyerType& destroyer,
        RepairerType& repairer,
        double time_limit_seconds,
        EvaluatorType evaluator = EvaluatorType{}) :
        m_destroyer(destroyer),
        m_repairer(repairer),
        m_evaluator(std::move(evaluator)),
        m_timeLimitSeconds(time_limit_seconds)
    {}

    template<typename SetCoverType, typename SetCoverContextType>
    auto operator()(std::shared_ptr<const SetCoverType> set_cover, SetCoverContextType context)
    {
        run(*set_cover, context);
        return std::tuple{std::move(set_cover), std::move(context)};
    }

    std::optional<StageMetrics> emit_metrics() const
    {
        if constexpr (RecordMetricsLevel > 0)
        {
            return m_metrics;
        }
        else
        {
            return std::nullopt;
        }
    }

    template<typename SetCoverType, typename SetCoverContextType>
    void run(const SetCoverType& set_cover, SetCoverContextType& context)
    {
        using Clock = std::chrono::steady_clock;
        const auto time_limit = std::chrono::duration<double>(m_timeLimitSeconds);
        const auto start = Clock::now();

        m_move.reserve(context.get_solution().size());
        while (Clock::now() - start < time_limit)
        {
            m_move.clear();
            SetCoverContextType candidate = context;
            m_destroyer.generateMove(set_cover, candidate, m_move);
            for (size_t set : m_move)
            {
                candidate.remove_set(set);
            }

            m_repairer.repair(set_cover, candidate);
            if (m_evaluator.evaluate(set_cover, context.get_solution(), candidate.get_solution()))
            {
                if constexpr (Debug)
                {
                    std::cout << "Improvement found! Cost: " << HeiConnect::sc::cost(set_cover, context.get_solution())
                              << " -> " << HeiConnect::sc::cost(set_cover, candidate.get_solution()) << std::endl;
                }
                context = std::move(candidate);
            }
        }

        if constexpr (RecordMetricsLevel > 0)
        {
            m_metrics = StageMetrics{
                {"cost", std::to_string(HeiConnect::sc::cost(set_cover, context.get_solution()))},
                {"size", std::to_string(context.get_solution().size())}
            };
        }
    }

private:
    DestroyerType& m_destroyer;
    RepairerType& m_repairer;
    EvaluatorType m_evaluator;
    double m_timeLimitSeconds;
    std::vector<size_t> m_move;
    
    std::optional<StageMetrics> m_metrics;
};

template<typename DestroyerType, typename RepairerType, size_t RecordMetricsLevel = 0>
BreakAndRepairSearch(DestroyerType&, RepairerType&, double, DefaultEvaluator = DefaultEvaluator{}) -> BreakAndRepairSearch<DestroyerType, RepairerType, RecordMetricsLevel>;


template<typename DestroyerType, typename RepairerType, typename EvaluatorType = DefaultEvaluator, bool Debug = false>
class LocalStepBreakAndRepairSearch
{
public:
    LocalStepBreakAndRepairSearch(
        RepairerType& repairer,
        double time_limit_seconds,
        EvaluatorType evaluator = EvaluatorType{},
        int max_steps = 10,
        std::mt19937_64 rng = std::mt19937_64{std::random_device{}()}) :
        m_repairer(repairer),
        m_evaluator(std::move(evaluator)),
        m_timeLimitSeconds(time_limit_seconds),
        m_nMaxSteps{max_steps},
        m_rng{rng}
    {}

    template<typename SetCoverType, typename ContextType>
    void run(const SetCoverType& set_cover, ContextType& context)
    {
        using Clock = std::chrono::steady_clock;
        const auto time_limit = std::chrono::duration<double>(m_timeLimitSeconds);
        const auto start = Clock::now();

        VectorSolution potential_sets{};
        std::vector<size_t> move_local;

        while (Clock::now() - start < time_limit)
        {
            move_local.clear();
            ContextType candidate = context;
            potential_sets.add_set(candidate.get_solution()[rand() % candidate.get_solution().size()]);
            for (int step{0}; step < m_nMaxSteps; step++)
            {
                size_t random_set = rand() % potential_sets.size();
                candidate.remove_set(potential_sets.get_solution()[random_set]);
                potential_sets.reset();
                BoundContext bound_context{candidate, potential_sets};
                m_repairer.repair(set_cover, bound_context);
            }

            if (m_evaluator.evaluate(set_cover, context.get_solution(), candidate.get_solution()))
            {
                if constexpr (Debug)
                {
                    std::cout << "Improvement found! Cost: " << HeiConnect::sc::cost(set_cover, context.get_solution())
                              << " -> " << HeiConnect::sc::cost(set_cover, candidate.get_solution()) << std::endl;
                }
                context = std::move(candidate);
            }
        }
    }

private:
    RepairerType& m_repairer;
    EvaluatorType m_evaluator;
    double m_timeLimitSeconds;
    int m_nMaxSteps;
    std::mt19937_64 m_rng;
};
