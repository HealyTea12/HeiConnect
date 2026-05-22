#pragma once
#include <chrono>
#include <vector>
#include <random>

#include "HeiConnect/set_cover/util.hpp"
#include "HeiConnect/set_cover/solver_greedy_context.hpp"


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

template<typename DestroyerType, typename RepairerType, typename EvaluatorType = DefaultEvaluator, bool Debug = false>
class BreakAndRepairSearch
{
public:
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

    template<typename SetCoverType, typename ContextType>
    void run(const SetCoverType& set_cover, ContextType& context)
    {
        using SetCost = typename SetCoverType::SetCost;
        using Clock = std::chrono::steady_clock;
        bool has_improvement = false;
        const auto time_limit = std::chrono::duration<double>(m_timeLimitSeconds);
        const auto start = Clock::now();

        m_move.reserve(context.get_solution().size());
        while (Clock::now() - start < time_limit)
        {
            m_move.clear();
            ContextType candidate = context;
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
    }

private:
    DestroyerType& m_destroyer;
    RepairerType& m_repairer;
    EvaluatorType m_evaluator;
    double m_timeLimitSeconds;
    std::vector<size_t> m_move;
};


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
        using SetCost = typename SetCoverType::SetCost;
        using Clock = std::chrono::steady_clock;
        bool has_improvement = false;
        const auto time_limit = std::chrono::duration<double>(m_timeLimitSeconds);
        const auto start = Clock::now();

        VectorSolution potential_sets{};

        std::vector<size_t> m_move_local;

        while (Clock::now() - start < time_limit)
        {
            m_move.clear();
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
    std::vector<size_t> m_move;
};