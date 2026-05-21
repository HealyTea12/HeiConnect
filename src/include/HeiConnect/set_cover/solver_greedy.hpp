#pragma once

#include "HeiConnect/set_cover/solver_base.hpp"
#include "HeiConnect/set_cover/solver_greedy_context.hpp"
#include "HeiConnect/set_cover/trimmer.hpp"

template <typename GreedyContext, typename SetCoverType>
concept GreedyContextCon = requires(GreedyContext context, size_t set_index) {
    { context.add_set(set_index) };
    { context.cover_count(set_index) } -> std::convertible_to<size_t>;
    { context.get_total_covered_elements() } -> std::convertible_to<size_t>;
};

// Generic greedy solver that works with any solver interface
template <typename SetCoverT, typename GreedyContext, typename Solution>
    requires GreedyContextCon<GreedyContext, SetCoverT>
class GreedyHelper
{
public:
    GreedyHelper() = default;

    void solve(const SetCoverT &set_cover, Solution &solution, GreedyContext &context)
    {
        std::priority_queue<std::pair<double, size_t>> pq;
        for (size_t s{0}; s < set_cover.get_num_sets(); s++)
        {
            size_t covered = 0;
            covered = context.cover_count(s);

            const double cost_benefit_ratio = static_cast<double>(covered) / set_cover.get_set_cost(s);
            pq.push({cost_benefit_ratio, s});
        }

        while (context.get_total_covered_elements() < set_cover.get_num_elements() &&
               solution.get_solution_size() != set_cover.get_num_sets())
        {
            auto best_set = pq.top();
            pq.pop();

            const auto best_set_idx = best_set.second;
            size_t covered = 0;
            covered = context.cover_count(best_set_idx);

            const double ratio = static_cast<double>(covered) / set_cover.get_set_cost(best_set_idx);
            if (ratio < best_set.first)
            {
                pq.push({ratio, best_set_idx});
                continue;
            }
            context.add_set(best_set_idx);
            solution.add_set(best_set_idx);
        }
    }
};

template <typename SetCoverType, typename SolutionType, typename GreedyContext = BasicContext<SetCoverType, SolutionType>>
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
        GreedyHelper<SetCoverType, GreedyContext, SolutionType> greedy{};
        greedy.solve(set_cover, solution, context);
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

template <typename SetCoverType, typename SolutionType>
class CheapestSetCoverSolver
{
public:
    void solve(const SetCoverType &set_cover, SolutionType &solution)
    {
        std::vector<size_t> set_indices(set_cover.get_num_sets());
        std::iota(set_indices.begin(), set_indices.end(), 0);
        std::sort(set_indices.begin(), set_indices.end(), [&](size_t a, size_t b)
        {
            return set_cover.get_set_cost(a) < set_cover.get_set_cost(b);
        });

        std::vector<size_t> coverage_counts(set_cover.get_num_elements(), 0);
        for (size_t set_index : set_indices)
        {
            if (std::all_of(coverage_counts.begin(), coverage_counts.end(), [](size_t count) { return count > 0; }))
            {
                break;
            }

            solution.add_set(set_index);
            set_cover.forEachElement(set_index, [&](size_t element)
            {
                ++coverage_counts[element];
            });
        }
    }
};
