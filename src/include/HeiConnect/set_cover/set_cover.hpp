#pragma once
#include <unordered_set>
#include <vector>
#include <stdexcept>

#include "HeiConnect/minmax.hpp"

struct SetCover
{
    const std::vector<size_t> a;
    const std::vector<size_t> b;
    const std::vector<double> costs;
    const size_t num_elements;

    SetCover(std::vector<size_t> a,
             std::vector<size_t> b,
             std::vector<double> costs,
             size_t num_elements)
        : a(std::move(a)), b(std::move(b)), costs(std::move(costs)), num_elements(num_elements)
    {
        if (this->a.size() - 1 != this->costs.size())
        {
            throw std::invalid_argument("SetCover: size of a and costs must be equal. Every subset should have a cost.");
        }
    };
    size_t get_num_elements() const noexcept
    {
        return num_elements;
    }
    double get_set_cost(size_t set_index) const
    {
        return costs[set_index];
    }
};

using ull = unsigned long long;
struct SetCoverBit
{
    const std::vector<ull> set_cover;
    const size_t n_sets;
    const size_t n_elements;
    const std::vector<double> costs;
    SetCoverBit(std::vector<ull> set_cover, size_t n_sets, size_t n_elements, std::vector<double> costs)
        : set_cover(std::move(set_cover)), n_sets(n_sets), n_elements(n_elements), costs(std::move(costs)) {};
    size_t get_num_elements() const noexcept
    {
        return n_elements;
    }
    double get_set_cost(size_t set_index) const
    {
        return costs[set_index];
    }
};

template <typename Derived, typename SetCoverType = SetCover>
    requires(std::same_as<SetCoverType, SetCover> || std::same_as<SetCoverType, SetCoverBit>)
class SetCoverSolver
{
public:
    SetCoverSolver(SetCoverType sc)
        : set_cover(std::move(sc))
    {
        NUM_ELEMENTS = set_cover.get_num_elements();
        m_covered_elements = std::vector<int>(NUM_ELEMENTS, 0);
    }
    ~SetCoverSolver() = default;
    void solve()
    {
        static_cast<Derived *>(this)->solve();
    };
    // Removes sets that are no longer necessary for the solution.
    // Those whose elements are still covered by other sets in the solution.
    void trim_solution()
    {
        while (true)
        {
            using WeightType = std::ranges::range_value_t<decltype(set_cover.costs)>;
            using SetType = std::ranges::range_value_t<decltype(set_cover.a)>;
            WeightType current_removable_weight{0};
            SetType current_removable_set{0};

            for (auto &set : m_solution)
            {
                bool can_remove = true;
                for (size_t j{set_cover.a[set]}; j < set_cover.a[set + 1]; j++)
                {
                    if (m_covered_elements[set_cover.b[j]] <= 1)
                    {
                        can_remove = false;
                        break;
                    }
                }
                if (can_remove > current_removable_weight)
                {
                    current_removable_weight = set_cover.costs[set];
                    current_removable_set = set;
                }
            }
            if (current_removable_weight > 0)
            {
                remove_set(current_removable_set);
            }
            else
            {
                break;
            }
        }
    }
    std::unordered_set<size_t> get_solution() const noexcept
    {
        return m_solution;
    }
    double get_solution_cost() const
    {
        double total_cost = 0.0;
        for (const auto &set_index : m_solution)
        {
            total_cost += set_cover.get_set_cost(set_index);
        }
        return total_cost;
    }

protected:
    void add_set(size_t set_index) noexcept
    {
        for (size_t j{set_cover.a[set_index]}; j < set_cover.a[set_index + 1]; j++)
        {
            if (m_covered_elements[set_cover.b[j]] == 0)
            {
                m_total_covered_elements += 1;
            }
            m_covered_elements[set_cover.b[j]] += 1;
        }
        m_solution.insert(set_index);
    }
    void remove_set(size_t set_index) noexcept
    {
        for (size_t j{set_cover.a[set_index]}; j < set_cover.a[set_index + 1]; j++)
        {
            m_covered_elements[set_cover.b[j]] -= 1;
            if (m_covered_elements[set_cover.b[j]] == 0)
            {
                m_total_covered_elements -= 1;
            }
        }
        m_solution.erase(set_index);
    }

protected:
    size_t NUM_ELEMENTS;
    SetCoverType set_cover;
    size_t m_total_covered_elements = 0;
    std::unordered_set<size_t> m_solution{};
    std::vector<bool> m_chosen_sets{};
    bool m_feasible = false;
    std::vector<int> m_covered_elements;
};

class SetCoverSolverGreedySingleThreadedPQ : public SetCoverSolver<SetCoverSolverGreedySingleThreadedPQ>
{
public:
    void solve();
};
class SetCoverSolverGreedyParallel : public SetCoverSolver<SetCoverSolverGreedyParallel>
{
public:
    void solve();
};

class SetCoverSolverSharpGreedy : public SetCoverSolver<SetCoverSolverSharpGreedy>
{
public:
    void solve();
};

class SetCoverSolverILP : public SetCoverSolver<SetCoverSolverILP>
{
public:
    using SetCoverSolver<SetCoverSolverILP>::SetCoverSolver;
    void solve();
};

class SetCoverSolverGreedyWideSingleThreaded : public SetCoverSolver<SetCoverSolverGreedyWideSingleThreaded>
{
public:
    void solve();
};