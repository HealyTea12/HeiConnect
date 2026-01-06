#pragma once
#include <vector>
#include <stdexcept>
#include <unordered_set>
#include <omp.h>
#include <bit>
#include <queue>
#include <algorithm>
#include <immintrin.h>

#include <gurobi_c++.h>

#include "HeiConnect/minmax.hpp"

template <typename T>
concept SetCoverCon = requires(T t, size_t i) {
    { t.get_num_sets() } -> std::convertible_to<size_t>;
    { t.get_num_elements() } -> std::convertible_to<size_t>;
    { t.get_set_cost(i) } -> std::floating_point;
    { t.set_begin(i) } -> std::input_iterator;
    { t.set_end(i) } -> std::input_iterator;
};

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
    size_t get_num_sets() const noexcept
    {
        return a.size() - 1;
    }
    size_t get_num_elements() const noexcept
    {
        return num_elements;
    }
    double get_set_cost(size_t set_index) const
    {
        return costs[set_index];
    }
    std::vector<size_t>::const_iterator set_begin(size_t set_index) const
    {
        return b.cbegin() + a[set_index];
    };
    std::vector<size_t>::const_iterator set_end(size_t set_index) const
    {
        return b.cbegin() + a[set_index + 1];
    }
};

static_assert(SetCoverCon<SetCover>);

using ull = unsigned long long;

class BitSetIterator
{
public:
    using value_type = size_t;
    using difference_type = std::ptrdiff_t;
    using iterator_concept = std::input_iterator_tag;
    using iterator_category = std::input_iterator_tag;

    BitSetIterator(const ull *data, size_t word_index, size_t n_words, bool is_end = false)
        : m_data(data), current_word_idx(word_index), current_word(0), m_is_end(is_end), m_n_words(n_words)
    {
        if (current_word_idx < m_n_words)
        {
            current_word = m_data[current_word_idx];
        }
        else
        {
            m_is_end = true;
        }
        next_bit();
    };
    BitSetIterator &operator++()
    {
        next_bit();
        return *this;
    }
    BitSetIterator operator++(int)
    {
        BitSetIterator tmp = *this;
        ++(*this);
        return tmp;
    }
    size_t operator*() const
    {
        return current_word_idx * sizeof(ull) * 8 + current_bit_idx;
    }
    bool operator!=(const BitSetIterator &other) const
    {
        // we only compare against the end
        return m_is_end != other.m_is_end;
        // return current_word_idx != other.current_word_idx || current_word != other.current_word;
    }
    friend bool operator==(const BitSetIterator &a, const BitSetIterator &b)
    {
        return !(a != b);
    }

private:
    void next_bit()
    {
        if (m_is_end)
            return;
        while (current_word == 0 && !m_is_end)
        {
            current_word_idx++;
            if (current_word_idx >= m_n_words)
            {
                m_is_end = true;
                return;
            }
            current_word = m_data[current_word_idx];
        }
        // find next set bit
        current_bit_idx = std::countr_zero(current_word);
        current_word ^= (1ULL << current_bit_idx);
    }
    size_t current_word_idx;
    ull current_word;
    size_t current_bit_idx = 0;
    const ull *m_data;
    bool m_is_end;
    size_t m_n_words;
};
static_assert(std::input_iterator<BitSetIterator>);

struct SetCoverBit
{
    const std::vector<ull> set_cover;
    const size_t n_sets;
    const size_t n_elements;
    const size_t n_cols;
    const std::vector<double> costs;
    SetCoverBit(std::vector<ull> set_cover, size_t n_sets, size_t n_elements, std::vector<double> costs)
        : set_cover(std::move(set_cover)), n_sets(n_sets), n_elements(n_elements), costs(std::move(costs)),
          n_cols((n_elements - 1) / (8 * sizeof(ull)) + 1) {};
    size_t get_num_sets() const noexcept
    {
        return n_sets;
    }
    size_t get_num_elements() const noexcept
    {
        return n_elements;
    }
    double get_set_cost(size_t set_index) const
    {
        return costs[set_index];
    }
    BitSetIterator set_begin(size_t set_index)
    {
        return BitSetIterator(set_cover.data() + set_index * n_cols, 0, n_cols);
    };
    BitSetIterator set_end(size_t set_index)
    {
        return BitSetIterator(set_cover.data() + set_index * n_cols, 0, n_cols, true); // only true matters
    }
};

static_assert(SetCoverCon<SetCoverBit>);

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
        using WeightType = decltype(set_cover.get_set_cost(0));
        using SetType = decltype(set_cover.get_num_sets()); // should be defined in SetCoverType
        while (true)
        {
            WeightType current_removable_weight{0};
            SetType current_removable_set{0};

            for (auto &set : m_solution)
            {
                bool can_remove = true;
                for (auto j{set_cover.set_begin(set)}; j != set_cover.set_end(set); j++)
                {
                    if (m_covered_elements[*j] == 1)
                    {
                        can_remove = false;
                        break;
                    }
                }
                if (can_remove && set_cover.get_set_cost(set) > current_removable_weight)
                {
                    current_removable_weight = set_cover.get_set_cost(set);
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
        for (auto j{set_cover.set_begin(set_index)}; j != set_cover.set_end(set_index); j++)
        {
            if (m_covered_elements[*j] == 0)
            {
                m_total_covered_elements += 1;
            }
            m_covered_elements[*j] += 1;
        }
        m_solution.insert(set_index);
    }
    void remove_set(size_t set_index) noexcept
    {
        for (auto j{set_cover.set_begin(set_index)}; j != set_cover.set_end(set_index); j++)
        {
            m_covered_elements[*j] -= 1;
            if (m_covered_elements[*j] == 0)
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

template <typename SetCoverType = SetCover>
class SetCoverSolverGreedyParallel : public SetCoverSolver<SetCoverSolverGreedyParallel<SetCoverType>, SetCoverType>
{
public:
    using Base = SetCoverSolver<SetCoverSolverGreedyParallel<SetCoverType>, SetCoverType>;
    using Base::add_set;
    using Base::m_covered_elements;
    using Base::m_solution;
    using Base::m_total_covered_elements;
    using Base::NUM_ELEMENTS;
    using Base::remove_set;
    using Base::set_cover;
    using Base::SetCoverSolver;

    void solve()
    {
        std::vector<double> cost_benefit_ratios = std::vector<double>(set_cover.a.size() - 1, 0.);
        // stop when all elements are covered or we've chosen all available sets
        while (m_total_covered_elements < NUM_ELEMENTS && m_solution.size() != set_cover.a.size() - 1)
        {
            // calculate ratios
#pragma omp parallel for
            for (size_t i = 0; i < set_cover.a.size() - 1; i++)
            {
                // skip already selected sets
                if (m_solution.find(i) != m_solution.end())
                {
                    cost_benefit_ratios[i] = -std::numeric_limits<double>::infinity();
                    continue;
                }
                int covered = 0;
                for (auto j = set_cover.set_begin(i); j != set_cover.set_end(i); j++)
                {
                    if (m_covered_elements[*j] == 0)
                        covered += 1;
                }
                if (covered > 0)
                    cost_benefit_ratios[i] = covered / set_cover.costs[i];
                else
                    cost_benefit_ratios[i] = 0;
            }
            auto best_set = argmax(cost_benefit_ratios);
            this->add_set(best_set);
        }
    };
};

// The idea of this method is to maintain a priority queue of sets based on their cost-benefit ration.
// Then, to select the next best set, we pop the queue. But we need to check if the ration has changed.
// The point is to recalculate the ratios as little as possible.
template <typename SetCoverType = SetCover>
class SetCoverSolverGreedySingleThreadedPQ : public SetCoverSolver<SetCoverSolverGreedySingleThreadedPQ<SetCoverType>, SetCoverType>
{
public:
    using Base = SetCoverSolver<SetCoverSolverGreedySingleThreadedPQ<SetCoverType>, SetCoverType>;
    using Base::add_set;
    using Base::m_covered_elements;
    using Base::m_solution;
    using Base::m_total_covered_elements;
    using Base::NUM_ELEMENTS;
    using Base::remove_set;
    using Base::set_cover;
    using Base::SetCoverSolver;

    void solve()
    {
        std::priority_queue<std::pair<double, size_t>> pq;
        // create a heap for the indices of all sets based on their cost-benefit ratio
        for (size_t i = 0; i < set_cover.get_num_sets(); i++)
        {
            int covered = 0;
            for (auto j = set_cover.set_begin(i); j != set_cover.set_end(i); j++)
            {
                if (m_covered_elements[*j] == 0)
                    covered += 1;
            }
            double ratio = (covered > 0) ? (covered / set_cover.costs[i]) : 0.;
            pq.push({ratio, i});
        }
        // stop when all elements are covered or we've chosen all available sets while (m_total_covered_elements < NUM_ELEMENTS && m_solution.size() != set_cover.a.size() - 1)
        while (m_total_covered_elements < NUM_ELEMENTS && m_solution.size() != set_cover.get_num_sets())
        {
            auto best_set = pq.top();
            pq.pop();
            // check if ratio has changed
            size_t best_set_idx = best_set.second;
            size_t covered = 0;
            for (auto i = set_cover.set_begin(best_set_idx); i != set_cover.set_end(best_set_idx); i++)
            {
                if (m_covered_elements[*i] == 0)
                    covered += 1;
            }
            double ratio = covered / set_cover.get_set_cost(best_set_idx);
            // if changed that put back in and try again
            if (ratio < best_set.first)
            {
                pq.push({ratio, best_set_idx});
                continue;
            }
            // otherwise add the set to the solution
            this->add_set(best_set_idx);
        }
    };
};

// Partial specialization for SetCoverBit
template <>
class SetCoverSolverGreedySingleThreadedPQ<SetCoverBit>
    : public SetCoverSolver<SetCoverSolverGreedySingleThreadedPQ<SetCoverBit>, SetCoverBit>
{
public:
    using Base = SetCoverSolver<SetCoverSolverGreedySingleThreadedPQ<SetCoverBit>, SetCoverBit>;
    using Base::add_set;
    using Base::m_covered_elements;
    using Base::m_solution;
    using Base::m_total_covered_elements;
    using Base::NUM_ELEMENTS;
    using Base::remove_set;
    using Base::set_cover;
    using Base::SetCoverSolver;

    size_t cover_count(const ull *set_data, const std::vector<ull> &covered_elements)
    {
        size_t covered = 0;
#ifdef __AVX512VPOPCNTDQ__
        for (size_t k{0}; k + 7 < set_cover.n_cols; k += 8)
        {
            __m512i vec_set = _mm512_loadu_si512(set_data + k);
            __m512i vec_covered = _mm512_loadu_si512(&covered_elements[k]);
            __m512i vec_new_bits = _mm512_andnot_si512(vec_covered, vec_set);
            __m512i vec_popcnt = _mm512_popcnt_epi64(vec_new_bits);
            covered += _mm512_reduce_add_epi64(vec_popcnt);
        }
        for (size_t k = (set_cover.n_cols / 8) * 8; k < set_cover.n_cols; k++)
        {
            ull new_bits = set_data[k] & ~covered_elements[k];
            covered += std::popcount(new_bits);
        }
#elif __AVX2__
        for (size_t k{0}; k + 3 < set_cover.n_cols; k += 4)
        {
            __m256i set_vec = _mm256_loadu_si256((__m256i *)(set_data + k));
            __m256i covered_vec = _mm256_loadu_si256((__m256i *)(&covered_elements[k]));
            __m256i new_bits = _mm256_andnot_si256(covered_vec, set_vec);
            covered += std::popcount(((ull *)&new_bits)[0]);
            covered += std::popcount(((ull *)&new_bits)[1]);
            covered += std::popcount(((ull *)&new_bits)[2]);
            covered += std::popcount(((ull *)&new_bits)[3]);
        }
        for (size_t k = (set_cover.n_cols / 4) * 4; k < set_cover.n_cols; k++)
        {
            ull new_bits = set_data[k] & ~covered_elements[k];
            covered += std::popcount(new_bits);
        }
#else
        for (size_t k{0}; k < set_cover.n_cols; k++)
        {
            ull new_bits = set_data[k] & ~covered_elements[k];
            covered += std::popcount(new_bits);
        }
#endif
        return covered;
    }

    void solve()
    {
        std::priority_queue<std::pair<double, size_t>> pq;
        std::vector<ull> covered_elements = std::vector<ull>(set_cover.n_cols, 0);
        covered_elements[set_cover.n_cols - 1] = 0xFFFFFFFFFFFFFFFFULL << (set_cover.get_num_elements() % (8 * sizeof(ull)));
        for (size_t i = 0; i < set_cover.get_num_sets(); i++)
        {
            const ull *set_data = &set_cover.set_cover[i * set_cover.n_cols];
            __builtin_prefetch(set_data, 0, 3);
            size_t covered = cover_count(set_data, covered_elements);
            double cost_benefit_ratio = static_cast<double>(covered) / set_cover.get_set_cost(i);
            pq.push({cost_benefit_ratio, i});
        }
        // stop when all elements are covered or we've chosen all available sets while (m_total_covered_elements < NUM_ELEMENTS && m_solution.size() != set_cover.a.size() - 1)
        while (m_total_covered_elements < set_cover.get_num_elements() && m_solution.size() != set_cover.get_num_sets())
        {
            auto best_set = pq.top();
            pq.pop();
            // check if ratio has changed
            size_t best_set_idx = best_set.second;
            const ull *set_data = &set_cover.set_cover[best_set_idx * set_cover.n_cols];
            __builtin_prefetch(set_data, 0, 3);
            size_t covered = cover_count(set_data, covered_elements);
            double ratio = static_cast<double>(covered) / set_cover.get_set_cost(best_set_idx);
            // if changed that put back in and try again
            if (ratio < best_set.first)
            {
                pq.push({ratio, best_set_idx});
                continue;
            }
            // otherwise add the set to the solution
            for (size_t k{0}; k < set_cover.n_cols; k++)
            {
                covered_elements[k] |= set_data[k];
            }
            m_total_covered_elements += covered;
            m_solution.insert(best_set_idx);
        }
    }
};

template <typename SetCoverType>
concept ElementIterableCon = requires(SetCoverType t, size_t e) {
    { t.element_begin(e) } -> std::input_iterator;
    { t.element_end(e) } -> std::input_iterator;
};

// assumes that a represents the elements, and b the sets covering them
// for this algorithm to work correctly, the represenation needs to be inverted, which is confusing
template <typename SetCoverType = SetCover>
    requires SetCoverCon<SetCoverType> && ElementIterableCon<SetCoverType>
class SetCoverSolverSharpGreedy : public SetCoverSolver<SetCoverSolverSharpGreedy<SetCoverType>, SetCoverType>
{
public:
    using Base = SetCoverSolver<SetCoverSolverSharpGreedy<SetCoverType>, SetCoverType>;
    using Base::add_set;
    using Base::m_covered_elements;
    using Base::m_solution;
    using Base::m_total_covered_elements;
    using Base::NUM_ELEMENTS;
    using Base::remove_set;
    using Base::set_cover;
    using Base::SetCoverSolver;

    void solve()
    {
        while (m_total_covered_elements < set_cover.get_num_elements() && m_solution.size() != set_cover.get_num_sets())
        {
            for (size_t e = 0; e < set_cover.get_num_elements(); e++)
            {
                if (m_covered_elements[e] > 0)
                    continue;
                auto cheapest_set_cost = std::numeric_limits<double>::max();
                size_t cheapest_set = std::numeric_limits<size_t>::max();
                // iterate over all sets containing element e (needs to be implemented)
                for (auto i = set_cover.element_begin(e); i != set_cover.element_end(e); i++)
                {
                    auto set = *i;
                    if (set_cover.get_set_cost(set) < cheapest_set_cost)
                    {
                        cheapest_set_cost = set_cover.get_set_cost(set);
                        cheapest_set = set;
                    }
                }
                if (cheapest_set == std::numeric_limits<size_t>::max())
                {
                    throw std::runtime_error(
                        "SetCoverSolverSharpGreedy: No set covers element " + std::to_string(e) + ". \n" +
                        "Problem instance is unsolvable.");
                }
                this->add_set(cheapest_set);
            }
        }
    }
};

template <typename SetCoverType = SetCover>
class SetCoverSolverGreedyCheapest : public SetCoverSolver<SetCoverSolverGreedyCheapest<SetCoverType>, SetCoverType>
{
public:
    using Base = SetCoverSolver<SetCoverSolverGreedyCheapest<SetCoverType>, SetCoverType>;
    using Base::add_set;
    using Base::m_chosen_sets;
    using Base::m_covered_elements;
    using Base::m_feasible;
    using Base::m_solution;
    using Base::m_total_covered_elements;
    using Base::NUM_ELEMENTS;
    using Base::remove_set;
    using Base::set_cover;
    using Base::SetCoverSolver;

    void solve()
    {
        while (m_total_covered_elements < set_cover.get_num_elements() && m_solution.size() < set_cover.get_num_sets())
        {
            double min_set_cost = std::numeric_limits<double>::max();
            size_t best_set = std::numeric_limits<size_t>::max();
            for (size_t s{0}; s < set_cover.get_num_sets(); s++)
            {
                if (m_solution.find(s) != m_solution.end())
                    continue;
                // check if it covers something
                bool covers_new = false;
                for (auto e = set_cover.set_begin(s); e != set_cover.set_end(s); e++)
                {
                    if (m_covered_elements[*e] == 0)
                    {
                        covers_new = true;
                        break;
                    }
                }
                if (!covers_new)
                    continue;
                double set_cost = set_cover.get_set_cost(s);
                if (set_cost < min_set_cost)
                {
                    min_set_cost = set_cost;
                    best_set = s;
                }
            }
            if (best_set == std::numeric_limits<size_t>::max())
            {
                m_feasible = false;
                return;
            }
            add_set(best_set);
        }
    }
};

template <typename SetCoverType = SetCover>
class SetCoverSolverILP : public SetCoverSolver<SetCoverSolverILP<SetCoverType>, SetCoverType>
{
public:
    using Base = SetCoverSolver<SetCoverSolverILP<SetCoverType>, SetCoverType>;
    using Base::add_set;
    using Base::m_chosen_sets;
    using Base::m_covered_elements;
    using Base::m_feasible;
    using Base::m_solution;
    using Base::m_total_covered_elements;
    using Base::NUM_ELEMENTS;
    using Base::remove_set;
    using Base::set_cover;
    using Base::SetCoverSolver;

    void solve()
    {
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "set_cover_ilp.log");
        env.start();
        GRBModel model = GRBModel(env);
        // create variables for each set
        std::vector<GRBVar> vars = std::vector<GRBVar>(set_cover.get_num_sets());
        for (size_t i = 0; i < set_cover.get_num_sets(); i++)
        {
            vars[i] = model.addVar(0.0, 1.0, set_cover.get_set_cost(i), GRB_BINARY, "s" + std::to_string(i));
        }
        // create coverage constraints: for each element e, sum_{sets i covering e} vars[i] >= 1
        // build linear expressions for each element by iterating sets
        std::vector<GRBLinExpr> cover_expr = std::vector<GRBLinExpr>(NUM_ELEMENTS);
        for (size_t i = 0; i < set_cover.get_num_sets(); i++)
        {
            for (auto j = set_cover.set_begin(i); j != set_cover.set_end(i); j++)
            {
                size_t elem = *j;
                cover_expr[elem] += vars[i];
            }
        }
        for (size_t e = 0; e < set_cover.get_num_elements(); e++)
        {
            model.addConstr(cover_expr[e] >= 1, "cover_e" + std::to_string(e));
        }
        // objective
        GRBLinExpr obj = 0;
        for (size_t i = 0; i < set_cover.get_num_sets(); i++)
        {
            obj += set_cover.get_set_cost(i) * vars[i];
        }
        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize();
        // check optimization status
        int status = model.get(GRB_IntAttr_Status);
        std::vector<bool> chosen_sets = std::vector<bool>(set_cover.get_num_sets(), false);
        if (status == GRB_OPTIMAL || status == GRB_SUBOPTIMAL)
        {
            for (size_t i = 0; i < set_cover.get_num_sets(); i++)
            {
                double val = vars[i].get(GRB_DoubleAttr_X);
                if (val > 0.5)
                    chosen_sets[i] = true;
            }
            m_feasible = true;
            m_chosen_sets = chosen_sets;
        }
        else
        {
            m_feasible = false;
        }
    }
};