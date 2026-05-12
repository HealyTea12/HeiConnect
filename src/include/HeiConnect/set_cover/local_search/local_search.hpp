#pragma once

#include <algorithm>
#include <cmath>
#include <omp.h>

template<typename MoveGeneratorType, typename RepairSolverType, double Eps = 1e-9, bool Debug = false>
class LocalSearchStage
{
public:
    LocalSearchStage(MoveGeneratorType& move_generator, RepairSolverType& repair_solver, double time_limit_seconds) :
        m_moveGenerator(move_generator),
        m_repairSolver(repair_solver),
        m_timeLimitSeconds(time_limit_seconds)
    {}

    template<typename SetCoverType, typename SolutionType, typename ContextType>
    void run(const SetCoverType& set_cover, SolutionType& solution, ContextType& context)
    {
        const double start_time = omp_get_wtime();

        double incumbent_cost = solution_cost(set_cover, solution);

        while (!time_exceeded(start_time))
        {
            bool improved = false;

            m_moveGenerator.forEachPotentialMove(set_cover, solution, [&](const auto& move) -> bool {
                if (time_exceeded(start_time))
                {
                    return false;
                }

                SolutionType candidate = solution;

                for (auto set : move)
                {
                    candidate.remove_set(set);
                }

                ContextType candidate_context = context;
                if constexpr (requires { candidate.remove_set(move[0]); })
                {
                    for (auto set : move)
                    {
                        candidate_context.remove_set(set);
                    }
                }
                else
                {
                    candidate_context.reset();
                    for (auto set : candidate.get_solution())
                    {
                        candidate_context.add_set(set);
                    }
                }


                repair(set_cover, candidate, candidate_context);

                const double candidate_cost = solution_cost(set_cover, candidate);

                if (candidate_cost + Eps < incumbent_cost)
                {
                    if constexpr (Debug)
                    {
                        std::cout << "Found improving move with cost " << candidate_cost
                                  << " (incumbent: " << incumbent_cost
                                  << ", improvement: " << incumbent_cost - candidate_cost << ")" << std::endl;
                    }
                    solution = std::move(candidate);
                    context = std::move(candidate_context);
                    incumbent_cost = candidate_cost;
                    improved = true;

                    // Stop enumeration and restart from the new solution.
                    return false;
                }

                return true;
            });

            if (!improved)
            {
                break;
            }
        }
    }

private:
    template<typename SetCoverType, typename SolutionType>
    static double solution_cost(const SetCoverType& set_cover, const SolutionType& solution)
    {
        double total = 0.0;

        for (auto set : solution.get_solution())
        {
            total += set_cover.get_set_cost(set);
        }

        return total;
    }

    bool time_exceeded(double start_time) const
    {
        return m_timeLimitSeconds > 0.0 && omp_get_wtime() - start_time >= m_timeLimitSeconds;
    }

    template<typename SetCoverType, typename SolutionType, typename ContextType>
    void repair(const SetCoverType& set_cover, SolutionType& candidate, ContextType& candidate_context)
    {
        if constexpr (requires { m_repairSolver.solve(set_cover, candidate, candidate_context); })
        {
            m_repairSolver.solve(set_cover, candidate, candidate_context);
        }
        else if constexpr (requires { m_repairSolver.run(set_cover, candidate, candidate_context); })
        {
            m_repairSolver.run(set_cover, candidate, candidate_context);
        }
        else if constexpr (requires { m_repairSolver.solve(set_cover, candidate); })
        {
            m_repairSolver.solve(set_cover, candidate);
        }
        else
        {
            m_repairSolver.run(set_cover, candidate);
        }
    }

    MoveGeneratorType& m_moveGenerator;
    RepairSolverType& m_repairSolver;
    double m_timeLimitSeconds;
};