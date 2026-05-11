#pragma once

#include "HeiConnect/set_cover/solver_base.hpp"
#include "HeiConnect/set_cover/solver_greedy_context.hpp"
#include "HeiConnect/set_cover/trimmer.hpp"

template <typename GreedyContext, typename SetCoverType>
concept GreedyContextCon = requires(GreedyContext context, size_t set_index, const SetCoverType &set_cover, USSolution solution) {
    { context.add_set(set_index) };
    { context.cover_count(set_index) } -> std::convertible_to<size_t>;
    { context.get_total_covered_elements() } -> std::convertible_to<size_t>;
    { set_cover.get_num_sets() } -> std::convertible_to<size_t>;
    { set_cover.get_set_cost(set_index) } -> std::convertible_to<double>;
    { set_cover.get_num_elements() } -> std::convertible_to<size_t>;
    { solution.add_set(set_index) };
    { solution.get_solution_size() } -> std::convertible_to<size_t>;
};

// Greedy solver with customizable context
template <typename SetCoverType, typename SolutionType, typename GreedyContext = BasicContext<SetCoverType, SolutionType>>
    requires GreedyContextCon<GreedyContext, SetCoverType>
class GreedySetCoverSolver
{
public:
    void solve(const SetCoverType &set_cover, SolutionType &solution)
    {
        auto set_cover_ptr = std::make_shared<SetCoverType>(set_cover);
        GreedyContext context{set_cover_ptr};
        solve(*set_cover_ptr, solution, context);
    }

    void solve(const SetCoverType &set_cover, SolutionType &solution, GreedyContext &context)
    {
        // Greedy algorithm: build priority queue of all sets by cost-benefit ratio
        std::priority_queue<std::pair<double, size_t>> pq;
        for (size_t s{0}; s < set_cover.get_num_sets(); s++)
        {
            const size_t covered = context.cover_count(s);
            const double cost_benefit_ratio = static_cast<double>(covered) / set_cover.get_set_cost(s);
            pq.push({cost_benefit_ratio, s});
        }

        // Greedily select sets until all elements are covered or all sets are selected
        while (context.get_total_covered_elements() < set_cover.get_num_elements() &&
               solution.get_solution_size() != set_cover.get_num_sets())
        {
            auto best_set = pq.top();
            pq.pop();

            const auto best_set_idx = best_set.second;
            const size_t covered = context.cover_count(best_set_idx);
            const double ratio = static_cast<double>(covered) / set_cover.get_set_cost(best_set_idx);

            // Re-evaluate if ratio has changed; if so, re-insert with new ratio
            if (ratio < best_set.first)
            {
                pq.push({ratio, best_set_idx});
                continue;
            }

            context.add_set(best_set_idx);
            solution.add_set(best_set_idx);
        }

        m_feasible = context.get_total_covered_elements() == set_cover.get_num_elements();
        m_total_covered_elements = context.get_total_covered_elements();
    }

    bool is_feasible() const noexcept
    {
        return m_feasible;
    }

    size_t get_total_covered_elements() const noexcept
    {
        return m_total_covered_elements;
    }

private:
    bool m_feasible = false;
    size_t m_total_covered_elements = 0;
};
