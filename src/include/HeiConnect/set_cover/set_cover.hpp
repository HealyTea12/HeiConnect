#pragma once
#include <vector>
#include <stdexcept>
#include <unordered_set>
#include <omp.h>
#include <bit>
#include <queue>
#include <algorithm>
#include <numeric>
#include <functional>
#include <immintrin.h>
#include <iomanip>
#include <cstdint>

#include <gurobi_c++.h>

#include "HeiConnect/minmax.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/set_cover/cactus_min_cuts.hpp"

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

template <class link_node_T, class link_edge_T>
    requires std::integral<link_node_T> && std::integral<link_edge_T>
class SetCoverPseudo
{
public:
    const std::vector<ull> min_cuts;
    const ull n_min_cuts;
    const std::vector<link_edge_T> link_vertices;
    const std::vector<link_node_T> link_edges;
    const std::vector<double> link_weights;
    const size_t n_cols;
    SetCoverPseudo(
        const std::vector<ull> min_cuts,
        ull n_min_cuts,
        const std::vector<link_edge_T> link_vertices,
        const std::vector<link_node_T> link_edges,
        const std::vector<double> link_weights)
        : min_cuts(std::move(min_cuts)),
          n_min_cuts(n_min_cuts),
          link_vertices(std::move(link_vertices)),
          link_edges(std::move(link_edges)),
          link_weights(std::move(link_weights)),
          n_cols(min_cuts.size() / (link_vertices.size() - 1)) {};
    size_t get_num_sets() const noexcept
    {
        return link_weights.size();
    }
    size_t get_num_elements() const noexcept
    {
        return n_min_cuts;
    }
    double get_set_cost(size_t set_index) const
    {
        return link_weights[set_index];
    }
    BitSetIterator set_begin(size_t set_index)
    {
        throw std::logic_error("SetCoverPseudo does not fully implement interface yet.");
    };
    BitSetIterator set_end(size_t set_index)
    {
        throw std::logic_error("SetCoverPseudo does not fully implement interface yet.");
    }
};

static_assert(SetCoverCon<SetCoverPseudo<uint64_t, uint64_t>>);

// TODO: set covers are assuming size_t indices everywhere, should improve with templates
template <class link_node_T, class link_edge_T, class link_weight_T>
    requires std::integral<link_node_T> && std::integral<link_edge_T>
class SetCoverOracle
{
public:
    SetCoverOracle(
        std::vector<uint32_t> tin,
        std::vector<uint32_t> tout,
        CactusMinCuts<size_t> mcs,
        std::vector<link_edge_T> link_vertices,
        std::vector<link_node_T> link_edges,
        std::vector<link_weight_T> link_weights) : m_tin(std::move(tin)), m_tout(std::move(tout)), m_mcs(std::move(mcs)),
                                                   link_vertices(std::move(link_vertices)), link_edges(std::move(link_edges)), link_weights(std::move(link_weights))
    {

        WeightedCRFGraph<> g = WeightedCRFGraph<>{
            {this->link_vertices, this->link_edges},
            this->link_weights};
        links = g.csr_to_vec_links();
    };

private:
    std::vector<uint32_t> m_tin;
    std::vector<uint32_t> m_tout;
    CactusMinCuts<size_t> m_mcs;

public:
    const std::vector<link_edge_T> link_vertices;
    const std::vector<link_node_T> link_edges;
    const std::vector<double> link_weights;
    std::vector<std::tuple<size_t, size_t, double>> links; // (u, v, weight)

public:
    size_t get_num_sets() const noexcept
    {
        return link_weights.size();
    }
    size_t get_num_elements() const noexcept
    {
        return m_mcs.get_n_min_cuts();
    }
    double get_set_cost(size_t set_index) const
    {
        return static_cast<double>(link_weights[set_index]);
    }
    BitSetIterator set_begin(size_t set_index)
    {
        throw std::logic_error("SetCoverPseudo does not fully implement interface yet.");
    };
    BitSetIterator set_end(size_t set_index)
    {
        throw std::logic_error("SetCoverPseudo does not fully implement interface yet.");
    }
    size_t get_num_tree_cuts() const
    {
        return m_mcs.TREE_CUTS.size();
    }
    size_t get_num_cycle_cuts() const
    {
        return m_mcs.CYCLE_CUTS.size();
    }
    bool covers_tree(size_t set_index, size_t element_index) const
    {
        auto [u, v, w] = links[set_index];
        auto cut = m_mcs.TREE_CUTS[element_index];
        bool is_ancestor_cut = isAncestor(cut, u) != isAncestor(cut, v);
        return is_ancestor_cut;
    }
    bool covers_cycle(size_t set_index, size_t element_index) const
    {
        auto [u, v, w] = links[set_index];
        auto cut = m_mcs.CYCLE_CUTS[element_index];
        bool u_belongs = isAncestor(cut.first, u) && !isAncestor(cut.second, u);
        bool v_belongs = isAncestor(cut.first, v) && !isAncestor(cut.second, v);
        return u_belongs != v_belongs;
    }

private:
    bool isAncestor(size_t anc, size_t node) const
    {
        return m_tin[anc] <= m_tin[node] && m_tout[anc] >= m_tout[node];
    }
};

template <typename Derived, typename SetCoverType = SetCover>
    requires(SetCoverCon<SetCoverType>)
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
        std::vector<size_t> solution_vec(m_solution.begin(), m_solution.end());
        // sort solution by decreasing cost
        std::sort(solution_vec.begin(), solution_vec.end(), [&](size_t a, size_t b)
                  { return set_cover.get_set_cost(a) > set_cover.get_set_cost(b); });
        for (const auto &set_index : solution_vec)
        {
            bool can_remove = true;
            for (auto j{set_cover.set_begin(set_index)}; j != set_cover.set_end(set_index); j++)
            {
                if (m_covered_elements[*j] == 1)
                {
                    can_remove = false;
                    break;
                }
            }
            if (can_remove)
            {
                remove_set(set_index);
            }
        }
    }

    /*
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
    */

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

// specialization of SetCoverSolver with SetCoverType = SetCoverBit
template <typename Derived>
class SetCoverSolver<Derived, SetCoverBit>
{
public:
    using SetCoverType = SetCoverBit;
    SetCoverSolver(SetCoverType sc)
        : set_cover(std::move(sc))
    {
        NUM_ELEMENTS = set_cover.get_num_elements();
        m_covered_elements = std::vector<ull>(set_cover.n_cols, 0);
        m_covered_elements[set_cover.n_cols - 1] = 0xFFFFFFFFFFFFFFFFULL << (set_cover.get_num_elements() % (8 * sizeof(ull)));
    }

    void solve()
    {
        static_cast<Derived *>(this)->solve();
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

    // must vectorize
    bool can_remove(size_t set_index)
    {
        for (size_t k{0}; k < set_cover.n_cols; k++)
        {
            ull set_coverage = set_cover.set_cover[set_index * set_cover.n_cols + k];

            for (const auto &other_set_index : m_solution)
            {
                if (other_set_index == set_index)
                    continue;
                set_coverage &= ~(set_cover.set_cover[other_set_index * set_cover.n_cols + k]);
                if (set_coverage == 0)
                    break;
            }
            if (set_coverage != 0)
            {
                return false;
            }
        }
        return true;
    }

    void init_coverage()
    {
        // set all covered elements to 0
        m_total_covered_elements = 0;
        m_covered_elements = std::vector<ull>(set_cover.n_cols, 0);
        m_covered_elements[set_cover.n_cols - 1] = 0xFFFFFFFFFFFFFFFFULL << (set_cover.get_num_elements() % (8 * sizeof(ull)));
    }

    void greedy_solve()
    {
        std::priority_queue<std::pair<double, size_t>> pq;
        for (size_t s{0}; s < set_cover.get_num_sets(); s++)
        {
            size_t covered = cover_count(s);
            double cost_benefit_ratio = static_cast<double>(covered) / set_cover.get_set_cost(s);
            pq.push({cost_benefit_ratio, s});
        }
        // stop when all elements are covered or we've chosen all available sets while (m_total_covered_elements < NUM_ELEMENTS && m_solution.size() != set_cover.a.size() - 1)
        while (m_total_covered_elements < set_cover.get_num_elements() && m_solution.size() != set_cover.get_num_sets())
        {
            auto best_set = pq.top();
            pq.pop();
            // check if ratio has changed
            auto best_set_idx = best_set.second;
            size_t covered = cover_count(best_set_idx);
            double ratio = static_cast<double>(covered) / set_cover.get_set_cost(best_set_idx);
            // if changed that put back in and try again
            if (ratio < best_set.first)
            {
                pq.push({ratio, best_set_idx});
                continue;
            }
            // otherwise add the set to the solution
            add_set(best_set_idx);
            m_total_covered_elements += covered;
            m_solution.insert(best_set_idx);
        }
    }

    // always the same could remove and just use base
    void trim_solution()
    {
        // sort solution by decreasing cost
        std::vector<size_t> solution_vec(m_solution.begin(), m_solution.end());
        std::sort(solution_vec.begin(), solution_vec.end(), [&](size_t a, size_t b)
                  { return set_cover.get_set_cost(a) > set_cover.get_set_cost(b); });
        for (const auto &set_index : solution_vec)
        {
            if (can_remove(set_index))
            {
                m_solution.erase(set_index);
            }
        }
    }

    void add_set(size_t set_index)
    {
        size_t k{0};
// check for 512 bit vectorization support and use it if available, otherwise check for 256 bit vectorization support and use it if available, otherwise do scalar
#if __AVX512__
        for (; k + 7 < set_cover.n_cols; k += 8)
        {
            __m512i set_vec = _mm512_loadu_si512((__m512i *)(set_cover.set_cover.data() + set_index * set_cover.n_cols + k));
            __m512i covered_vec = _mm512_loadu_si512((__m512i *)(&m_covered_elements[k]));
            __m512i new_coverage = _mm512_or_si512(set_vec, covered_vec);
            _mm512_storeu_si512((__m512i *)(&m_covered_elements[k]), new_coverage);
        }
#elif __AVX2__
        for (; k + 3 < set_cover.n_cols; k += 4)
        {
            __m256i set_vec = _mm256_loadu_si256((__m256i *)(set_cover.set_cover.data() + set_index * set_cover.n_cols + k));
            __m256i covered_vec = _mm256_loadu_si256((__m256i *)(&m_covered_elements[k]));
            __m256i new_coverage = _mm256_or_si256(set_vec, covered_vec);
            _mm256_storeu_si256((__m256i *)(&m_covered_elements[k]), new_coverage);
        }
#endif
        for (; k < set_cover.n_cols; k++)
        {
            m_covered_elements[k] |= set_cover.set_cover[set_index * set_cover.n_cols + k];
        }
    }

    // must check instructions
    bool check_solved(std::vector<ull> &covered_els)
    {
        size_t k{0};
#if __AVX512__
        for (; k + 7 < set_cover.n_cols; k += 8)
        {
            __m512i covered_vec = _mm512_loadu_si512((__m512i *)(&covered_els[k]));
            __m512i all_covered = _mm512_set1_epi64(0xFFFFFFFFFFFFFFFFULL);
            __m512i cmp = _mm512_cmpeq_epi64_mask(covered_vec, all_covered);
            if (cmp != 0xFF)
            {
                return false;
            }
        }
#elif __AVX2__
        for (; k + 3 < set_cover.n_cols; k += 4)
        {
            __m256i covered_vec = _mm256_loadu_si256((__m256i *)(&covered_els[k]));
            __m256i all_covered = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFFULL);
            __m256i cmp = _mm256_cmpeq_epi64(covered_vec, all_covered);
            int mask = _mm256_movemask_pd(_mm256_castsi256_pd(cmp));
            if (mask != 0xF)
            {
                return false;
            }
        }
#endif
        for (; k < set_cover.n_cols; k++)
        {
            if (covered_els[k] != 0xFFFFFFFFFFFFFFFFULL)
            {
                return false;
            }
        }
        return true;
    }

    size_t cover_count(size_t set_index)
    {
        size_t n_covered_els = 0;
        size_t k{0};
#ifdef __AVX512VPOPCNTDQ__
        for (; k + 7 < set_cover.n_cols; k += 8)
        {
            __m512i vec_set = _mm512_loadu_si512((__m512i *)(set_cover.set_cover.data() + set_index * set_cover.n_cols + k));
            __m512i vec_covered = _mm512_loadu_si512(&m_covered_elements[k]);
            __m512i vec_new_bits = _mm512_andnot_si512(vec_covered, vec_set);
            __m512i vec_popcnt = _mm512_popcnt_epi64(vec_new_bits);
            n_covered_els += _mm512_reduce_add_epi64(vec_popcnt);
        }
#elif __AVX512__
        for (; k + 7 < set_cover.n_cols; k += 8)
        {
            __m512i vec_set = _mm512_loadu_si512((__m512i *)(set_cover.set_cover.data() + set_index * set_cover.n_cols + k));
            __m512i vec_covered = _mm512_loadu_si512(&m_covered_elements[k]);
            __m512i vec_new_bits = _mm512_andnot_si512(vec_covered, vec_set);
            // there is no popcnt instruction in AVX512 without VPOPCNTDQ, so we have to do it manually
            for (int i = 0; i < 8; i++)
            {
                n_covered_els += std::popcount(((ull *)&vec_new_bits)[i]);
            }
        }
#elif __AVX2__
        for (; k + 3 < set_cover.n_cols; k += 4)
        {
            __m256i vec_set = _mm256_loadu_si256((__m256i *)(set_cover.set_cover.data() + set_index * set_cover.n_cols + k));
            __m256i covered_vec = _mm256_loadu_si256((__m256i *)(&m_covered_elements[k]));
            __m256i new_bits = _mm256_andnot_si256(covered_vec, vec_set);
            for (int i = 0; i < 4; i++)
            {
                n_covered_els += std::popcount(((ull *)&new_bits)[i]);
            }
        }
#endif
        {
            const ull *set_data = set_cover.set_cover.data() + set_index * set_cover.n_cols;
            for (; k < set_cover.n_cols; k++)
            {
                ull covered_cuts = set_data[k] & ~m_covered_elements[k];
                n_covered_els += std::popcount(covered_cuts);
            }
        }
        return n_covered_els;
    }

protected:
    size_t NUM_ELEMENTS;
    SetCoverType set_cover;
    size_t m_total_covered_elements = 0;
    std::unordered_set<size_t> m_solution{};
    std::vector<ull> m_covered_elements;
    // ILP stuff don't really belong here
    std::vector<bool> m_chosen_sets{};
    bool m_feasible = false;
};

// specialization of SetCoverSolver  with SetCoverType = SetCoverPseudo
template <typename Derived, class link_node_T, class link_edge_T>
class SetCoverSolver<Derived, SetCoverPseudo<link_node_T, link_edge_T>>
{
public:
    using SetCoverType = SetCoverPseudo<link_node_T, link_edge_T>;

    SetCoverSolver(SetCoverType sc)
        : set_cover(std::move(sc))
    {
        NUM_ELEMENTS = set_cover.get_num_elements();
        m_covered_elements = std::vector<ull>(set_cover.n_cols, 0);
        m_covered_elements[set_cover.n_cols - 1] = 0xFFFFFFFFFFFFFFFFULL << (set_cover.get_num_elements() % (8 * sizeof(ull)));
        WeightedCRFGraph<> g = WeightedCRFGraph<>{
            {set_cover.link_vertices, set_cover.link_edges},
            set_cover.link_weights};
        m_links = g.csr_to_vec_links();
    }

    ~SetCoverSolver() = default;

    void solve()
    {
        static_cast<Derived *>(this)->solve();
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

    bool can_remove(size_t set_index)
    {
        auto [u, v, w] = m_links[set_index];
        for (size_t k{0}; k < set_cover.n_cols; k++)
        {
            ull set_coverage = set_cover.min_cuts[u * set_cover.n_cols + k] ^
                               set_cover.min_cuts[v * set_cover.n_cols + k];
            for (const auto &other_set_index : m_solution)
            {
                if (other_set_index == set_index)
                    continue;
                auto [ou, ov, ow] = m_links[other_set_index];
                set_coverage &= ~(set_cover.min_cuts[ou * set_cover.n_cols + k] ^
                                  set_cover.min_cuts[ov * set_cover.n_cols + k]);
                if (set_coverage == 0)
                    break;
            }
            if (set_coverage != 0)
            {
                return false;
            }
        }
        return true;
    }

    void init_coverage()
    {
        // set all covered elements to 0
        m_total_covered_elements = 0;
        m_covered_elements = std::vector<ull>(set_cover.n_cols, 0);
        m_covered_elements[set_cover.n_cols - 1] = 0xFFFFFFFFFFFFFFFFULL << (set_cover.get_num_elements() % (8 * sizeof(ull)));
    }

    void greedy_solve()
    {
        std::priority_queue<std::pair<double, std::tuple<size_t, size_t, size_t>>> pq;
        for (size_t u{0}; u < set_cover.link_vertices.size() - 1; u++)
        {
            for (size_t e{set_cover.link_vertices[u]}; e < set_cover.link_vertices[u + 1]; e++)
            {
                ull v = set_cover.link_edges[e];
                size_t covered = cover_count(u, v, m_covered_elements);
                double cost_benefit_ratio = static_cast<double>(covered) / set_cover.get_set_cost(e);
                pq.push({cost_benefit_ratio, {u, v, e}});
            }
        }
        // stop when all elements are covered or we've chosen all available sets while (m_total_covered_elements < NUM_ELEMENTS && m_solution.size() != set_cover.a.size() - 1)
        while (m_total_covered_elements < set_cover.get_num_elements() && m_solution.size() != set_cover.get_num_sets())
        {
            auto best_set = pq.top();
            pq.pop();
            // check if ratio has changed
            auto best_set_idx = best_set.second;
            size_t covered = cover_count(std::get<0>(best_set_idx), std::get<1>(best_set_idx), m_covered_elements);
            double ratio = static_cast<double>(covered) / set_cover.get_set_cost(std::get<2>(best_set_idx));
            // if changed that put back in and try again
            if (ratio < best_set.first)
            {
                pq.push({ratio, best_set_idx});
                continue;
            }
            // otherwise add the set to the solution
            add_set(m_covered_elements, std::get<0>(best_set_idx), std::get<1>(best_set_idx));
            m_total_covered_elements += covered;
            m_solution.insert(std::get<2>(best_set_idx));
        }
    }

    void trim_solution()
    {
        // sort solution by decreasing cost
        std::vector<size_t> solution_vec(m_solution.begin(), m_solution.end());
        std::sort(solution_vec.begin(), solution_vec.end(), [&](size_t a, size_t b)
                  { return set_cover.get_set_cost(a) > set_cover.get_set_cost(b); });
        for (const auto &set_index : solution_vec)
        {
            if (can_remove(set_index))
            {
                m_solution.erase(set_index);
            }
        }
    }

    // TODO: Doesn't belong here and could return an iterator
    // creates a vector of all subsets of size k from n elements, used for local search
    // for example k=2, n=4 -> {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}}
    // so returns a vector  -> [0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3]
    std::vector<size_t> calculate_subsets(size_t k, size_t n)
    {
        std::vector<size_t> subsets;
        std::vector<size_t> current_subset(k);
        std::function<void(size_t, size_t)> backtrack = [&](size_t start, size_t depth)
        {
            if (depth == k)
            {
                subsets.insert(subsets.end(), current_subset.begin(), current_subset.end());
                return;
            }
            for (size_t i = start; i < n; i++)
            {
                current_subset[depth] = i;
                backtrack(i + 1, depth + 1);
            }
        };
        backtrack(0, 0);
        return subsets;
    }

    // removes k sets from the current solution and tries to re-add better ones
    bool local_search(int k)
    {
        auto current_solution_cost = get_solution_cost();
        std::vector<size_t> solution_vec = std::vector<size_t>(m_solution.begin(), m_solution.end());
        // maybe this should be an iterator, because storing all this crap is unnecessary and inconvenient
        auto subsets = calculate_subsets(k, solution_vec.size());
        auto original_solution{m_solution};
        for (size_t i{0}; i < subsets.size(); i += k)
        {
            // save original solution

            // remove k sets
            for (size_t j{0}; j < k; j++)
            {
                size_t set_index = solution_vec[subsets[i + j]];
                m_solution.erase(set_index);
            }
            // refill covered elements
            init_coverage();
            for (const auto &other_set_index : m_solution)
            {
                auto [ou, ov, ow] = m_links[other_set_index];
                add_set(m_covered_elements, ou, ov);
            }
            m_total_covered_elements = calculate_covered_elements();
            // try to add k better sets
            greedy_solve();
            auto new_solution_cost = get_solution_cost();
            constexpr double epsilon = 1e-9;
            if (new_solution_cost < current_solution_cost - epsilon)
            {
                std::cout << k << " local search:" << "Improved cost from "
                          << std::fixed << std::setprecision(10) << current_solution_cost << " to " << new_solution_cost << std::endl;
                return true;
            }
            else
            {
                // revert to original solution (removes sets added by greedy_solve)
                m_solution = original_solution;
            }
        }
        return false;
    }

    size_t calculate_covered_elements()
    {
        size_t total_covered = 0;
        for (size_t k{0}; k < set_cover.n_cols; k++)
        {
            ull col = m_covered_elements[k];
            total_covered += std::popcount(col);
        }
        size_t remainder = set_cover.get_num_elements() % (8 * sizeof(ull));
        size_t tail = remainder ? (8 * sizeof(ull) - remainder) : 0;
        return total_covered - tail;
    }

    void add_set(std::vector<ull> &covered_els, ull u, ull v)
    {
        for (size_t k{0}; k < set_cover.n_cols; k++)
        {
            covered_els[k] |= set_cover.min_cuts[u * set_cover.n_cols + k] ^
                              set_cover.min_cuts[v * set_cover.n_cols + k];
        }
    }

    bool check_solved(std::vector<ull> &covered_els)
    {
        for (size_t k{0}; k < set_cover.n_cols; k++)
        {
            if (covered_els[k] != 0xFFFFFFFFFFFFFFFFULL)
            {
                return false;
            }
        }
        return true;
    }

    size_t cover_count(ull u, ull v, const std::vector<ull> &covered_elements)
    {
        const ull *u_data = &set_cover.min_cuts[u * set_cover.n_cols];
        const ull *v_data = &set_cover.min_cuts[v * set_cover.n_cols];
        size_t n_covered_cuts = 0;
#ifdef __AVX512VPOPCNTDQ__
        for (size_t k{0}; k + 7 < set_cover.n_cols; k += 8)
        {
            __m512i vec_u = _mm512_loadu_si512(u_data + k);
            __m512i vec_v = _mm512_loadu_si512(v_data + k);
            __m512i vec_covered = _mm512_loadu_si512(&covered_elements[k]);
            __m512i vec_covered_cuts = _mm512_xor_si512(vec_u, vec_v);
            __m512i vec_new_bits = _mm512_andnot_si512(vec_covered, vec_covered_cuts);
            __m512i vec_popcnt = _mm512_popcnt_epi64(vec_new_bits);
            n_covered_cuts += _mm512_reduce_add_epi64(vec_popcnt);
        }
        for (size_t k = (set_cover.n_cols / 8) * 8; k < set_cover.n_cols; k++)
        {
            ull covered_cuts = u_data[k] ^ v_data[k];
            n_covered_cuts += std::popcount(covered_cuts & ~covered_elements[k]);
        }
#elif __AVX2__
        for (size_t k{0}; k + 3 < set_cover.n_cols; k += 4)
        {
            __m256i u_vec = _mm256_loadu_si256((__m256i *)(u_data + k));
            __m256i v_vec = _mm256_loadu_si256((__m256i *)(v_data + k));
            __m256i covered_vec = _mm256_loadu_si256((__m256i *)(&covered_elements[k]));
            __m256i covered_cuts = _mm256_xor_si256(u_vec, v_vec);
            __m256i new_bits = _mm256_andnot_si256(covered_vec, covered_cuts);
            n_covered_cuts += std::popcount(((ull *)&new_bits)[0]);
            n_covered_cuts += std::popcount(((ull *)&new_bits)[1]);
            n_covered_cuts += std::popcount(((ull *)&new_bits)[2]);
            n_covered_cuts += std::popcount(((ull *)&new_bits)[3]);
        }
        for (size_t k = (set_cover.n_cols / 4) * 4; k < set_cover.n_cols; k++)
        {
            ull covered_cuts = u_data[k] ^ v_data[k];
            n_covered_cuts += std::popcount(covered_cuts & ~covered_elements[k]);
        }
#else
        for (size_t k{0}; k < set_cover.n_cols; k++)
        {
            ull covered_cuts = u_data[k] ^ v_data[k];
            n_covered_cuts += std::popcount(covered_cuts & ~covered_elements[k]);
        }
#endif
        return n_covered_cuts;
    }

protected:
    size_t NUM_ELEMENTS;
    SetCoverType set_cover;
    size_t m_total_covered_elements = 0;
    std::unordered_set<size_t> m_solution{};
    std::vector<bool> m_chosen_sets{};
    bool m_feasible = false;
    std::vector<ull> m_covered_elements;
    std::vector<std::tuple<size_t, size_t, double>> m_links;
};

// specialization of SetCoverSolver  with SetCoverType = SetCoverOracle
template <typename Derived, class link_node_T, class link_edge_T, class link_weight_T>
class SetCoverSolver<Derived, SetCoverOracle<link_node_T, link_edge_T, link_weight_T>>
{
public:
    using SetCoverType = SetCoverOracle<link_node_T, link_edge_T, link_weight_T>;

    SetCoverSolver(SetCoverType sc)
        : set_cover(std::move(sc))
    {
        NUM_ELEMENTS = set_cover.get_num_elements();
        m_coverage_count = std::vector<ull>(set_cover.get_num_elements(), 0);
        WeightedCRFGraph<> g = WeightedCRFGraph<>{
            {set_cover.link_vertices, set_cover.link_edges},
            set_cover.link_weights};
        m_links = g.csr_to_vec_links();
    }

    ~SetCoverSolver() = default;

    void solve()
    {
        static_cast<Derived *>(this)->solve();
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

    bool can_remove(size_t set_index)
    {
        auto [u, v, w] = m_links[set_index];
        for (size_t k{0}; k < set_cover.get_num_tree_cuts(); k++)
        {
            if (set_cover.covers_tree(set_index, k) && m_coverage_count[k] == 1)
            {
                return false;
            }
        }
        for (size_t k{0}; k < set_cover.get_num_cycle_cuts(); k++)
        {
            if (set_cover.covers_cycle(set_index, k) && m_coverage_count[set_cover.get_num_tree_cuts() + k] == 1)
            {
                return false;
            }
        }
        return true;
    }

    void greedy_solve()
    {
        std::priority_queue<std::pair<double, size_t>> pq;
        for (size_t s{0}; s < set_cover.get_num_sets(); s++)
        {
            size_t covered = cover_count(s);
            double cost_benefit_ratio = static_cast<double>(covered) / set_cover.get_set_cost(s);
            pq.push({cost_benefit_ratio, s});
        }
        // stop when all elements are covered or we've chosen all available sets while (m_total_covered_elements < NUM_ELEMENTS && m_solution.size() != set_cover.a.size() - 1)
        while (m_total_covered_elements < set_cover.get_num_elements() && m_solution.size() != set_cover.get_num_sets())
        {
            auto best_set = pq.top();
            pq.pop();
            // check if ratio has changed
            auto best_set_idx = best_set.second;
            size_t covered = cover_count(best_set_idx);
            double ratio = static_cast<double>(covered) / set_cover.get_set_cost(best_set_idx);
            // if changed that put back in and try again
            if (ratio < best_set.first)
            {
                pq.push({ratio, best_set_idx});
                continue;
            }
            // otherwise add the set to the solution
            add_set(best_set_idx);
        }
    }

    void trim_solution()
    {
        // sort solution by decreasing cost
        std::vector<size_t> solution_vec(m_solution.begin(), m_solution.end());
        std::sort(solution_vec.begin(), solution_vec.end(), [&](size_t a, size_t b)
                  { return set_cover.get_set_cost(a) > set_cover.get_set_cost(b); });
        for (const auto &set_index : solution_vec)
        {
            if (can_remove(set_index))
            {
                delete_set(set_index);
            }
        }
    }

    size_t calculate_covered_elements()
    {
        return std::count_if(m_coverage_count.begin(), m_coverage_count.end(), [](ull count)
                             { return count > 0; });
    }

    size_t cover_count(size_t set_index)
    {
        size_t covered = 0;
        auto [u, v, w] = m_links[set_index];
        for (size_t k{0}; k < set_cover.get_num_tree_cuts(); k++)
        {
            if (set_cover.covers_tree(set_index, k) && m_coverage_count[k] == 0)
            {
                covered++;
            }
        }
        for (size_t k{0}; k < set_cover.get_num_cycle_cuts(); k++)
        {
            if (set_cover.covers_cycle(set_index, k) && m_coverage_count[set_cover.get_num_tree_cuts() + k] == 0)
            {
                covered++;
            }
        }
        return covered;
    }

    void add_set(size_t set_index)
    {
        auto [u, v, w] = m_links[set_index];
        for (size_t k{0}; k < set_cover.get_num_tree_cuts(); k++)
        {
            if (set_cover.covers_tree(set_index, k))
            {
                size_t cov = m_coverage_count[k]++;
                if (cov == 0)
                    m_total_covered_elements++;
            }
        }
        for (size_t k{0}; k < set_cover.get_num_cycle_cuts(); k++)
        {
            if (set_cover.covers_cycle(set_index, k))
            {
                m_coverage_count[set_cover.get_num_tree_cuts() + k]++;
                if (m_coverage_count[set_cover.get_num_tree_cuts() + k] == 1)
                    m_total_covered_elements++;
            }
        }
        m_solution.insert(set_index);
    }

    void delete_set(size_t set_index)
    {
        auto [u, v, w] = m_links[set_index];
        for (size_t k{0}; k < set_cover.get_num_tree_cuts(); k++)
        {
            if (set_cover.covers_tree(set_index, k))
            {
                size_t cov = --m_coverage_count[k];
                if (cov == 0)
                    m_total_covered_elements--;
            }
        }
        for (size_t k{0}; k < set_cover.get_num_cycle_cuts(); k++)
        {
            if (set_cover.covers_cycle(set_index, k))
            {
                size_t cov = --m_coverage_count[set_cover.get_num_tree_cuts() + k];
                if (cov == 0)
                    m_total_covered_elements--;
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
    std::vector<ull> m_covered_elements;
    std::vector<ull> m_coverage_count;
    std::vector<std::tuple<size_t, size_t, double>> m_links;
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
            double ratio = (covered > 0) ? (covered / set_cover.get_set_cost(i)) : 0.;
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
    using Base::set_cover;
    using Base::SetCoverSolver;

    void solve()
    {
        init_coverage();
        greedy_solve();
    }
};

// Partial specialization for SetCoverPseudo
template <class link_node_T, class link_edge_T>
class SetCoverSolverGreedySingleThreadedPQ<SetCoverPseudo<link_node_T, link_edge_T>>
    : public SetCoverSolver<SetCoverSolverGreedySingleThreadedPQ<SetCoverPseudo<link_node_T, link_edge_T>>, SetCoverPseudo<link_node_T, link_edge_T>>
{
public:
    using Base = SetCoverSolver<SetCoverSolverGreedySingleThreadedPQ<SetCoverPseudo<link_node_T, link_edge_T>>, SetCoverPseudo<link_node_T, link_edge_T>>;
    using Base::add_set;
    using Base::m_covered_elements;
    using Base::m_total_covered_elements;
    using Base::NUM_ELEMENTS;
    using Base::set_cover;
    using Base::SetCoverSolver;
    // std::unordered_set<std::pair<size_t, size_t>> m_solution{};

    // we should reuse this code and avoid copying it
    void solve()
    {
        this->init_coverage();
        this->greedy_solve();
    }
};

// Partial specialization for SetCoverOracle
template <class link_node_T, class link_edge_T, class link_weight_T>
class SetCoverSolverGreedySingleThreadedPQ<SetCoverOracle<link_node_T, link_edge_T, link_weight_T>>
    : public SetCoverSolver<SetCoverSolverGreedySingleThreadedPQ<SetCoverOracle<link_node_T, link_edge_T, link_weight_T>>, SetCoverOracle<link_node_T, link_edge_T, link_weight_T>>
{
public:
    using Base = SetCoverSolver<SetCoverSolverGreedySingleThreadedPQ<SetCoverOracle<link_node_T, link_edge_T, link_weight_T>>, SetCoverOracle<link_node_T, link_edge_T, link_weight_T>>;

    // we should reuse this code and avoid copying it
    void solve()
    {
        Base::greedy_solve();
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
        std::vector<size_t> set_indices = std::vector<size_t>(set_cover.get_num_sets(), 0);
        std::iota(set_indices.begin(), set_indices.end(), 0);
        std::sort(set_indices.begin(), set_indices.end(), [&](size_t a, size_t b)
                  { return set_cover.get_set_cost(a) < set_cover.get_set_cost(b); });
        size_t i = 0;
        while (m_total_covered_elements < set_cover.get_num_elements() && m_solution.size() < set_cover.get_num_sets())
        {
            add_set(set_indices[i]);
            i++;
        }
    }

    /*
        void solve()
        {
            while (m_total_covered_elements < set_cover.get_num_elements() && m_solution.size() < set_cover.get_num_sets())
            {
                double min_set_cost = std::numeric_limits<double>::max();
                size_t best_set = std::numeric_limits<size_t>::max();
                for (size_t s{0}; s < set_cover.get_num_sets(); s++)
                {
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
    */
};

// Partial specialization for GreedyCheapest with SetCoverBit
template <>
class SetCoverSolverGreedyCheapest<SetCoverBit>
    : public SetCoverSolver<SetCoverSolverGreedyCheapest<SetCoverBit>, SetCoverBit>
{
public:
    using Base = SetCoverSolver<SetCoverSolverGreedyCheapest<SetCoverBit>, SetCoverBit>;
    using Base::add_set;
    using Base::m_chosen_sets;
    using Base::m_feasible;
    using Base::m_solution;
    using Base::m_total_covered_elements;
    using Base::NUM_ELEMENTS;
    using Base::set_cover;
    using Base::SetCoverSolver;

    using Base::m_covered_elements;

    // TODO: this is the same as for pseudo
    void solve()
    {
        init_coverage();
        std::vector<size_t> set_indices = std::vector<size_t>(set_cover.get_num_sets());
        std::iota(set_indices.begin(), set_indices.end(), 0);
        std::sort(set_indices.begin(), set_indices.end(), [&](size_t a, size_t b)
                  { return set_cover.get_set_cost(a) < set_cover.get_set_cost(b); });
        bool solved = false;
        size_t i = 0;
        while (!solved && m_solution.size() < set_cover.get_num_sets())
        {
            add_set(set_indices[i]);
            solved = check_solved(m_covered_elements);
            m_solution.insert(set_indices[i]);
            i++;
        }
    }
};

// Partial specialization for GreedyCheapest with SetCoverPseudo
template <typename link_node_T, typename link_edge_T>
class SetCoverSolverGreedyCheapest<SetCoverPseudo<link_node_T, link_edge_T>>
    : public SetCoverSolver<SetCoverSolverGreedyCheapest<SetCoverPseudo<link_node_T, link_edge_T>>, SetCoverPseudo<link_node_T, link_edge_T>>
{
public:
    using Base = SetCoverSolver<SetCoverSolverGreedyCheapest<SetCoverPseudo<link_node_T, link_edge_T>>, SetCoverPseudo<link_node_T, link_edge_T>>;
    using Base::add_set;
    using Base::m_chosen_sets;
    using Base::m_feasible;
    using Base::m_solution;
    using Base::m_total_covered_elements;
    using Base::NUM_ELEMENTS;
    using Base::set_cover;
    using Base::SetCoverSolver;

    using Base::m_covered_elements;
    using Base::m_links;

    void solve()
    {
        double start = omp_get_wtime();
        auto csr = WeightedCRFGraph(set_cover.link_vertices,
                                    set_cover.link_edges,
                                    set_cover.link_weights);
        m_links = csr.csr_to_vec_links();
        double end = omp_get_wtime();
        std::cout << "Converted CSR to vector of links in " << end - start << " seconds." << std::endl;
        m_covered_elements = std::vector<ull>(set_cover.n_cols, 0);
        m_covered_elements[set_cover.n_cols - 1] = 0xFFFFFFFFFFFFFFFFULL << (set_cover.get_num_elements() % (8 * sizeof(ull)));
        std::vector<size_t> set_indices = std::vector<size_t>(set_cover.get_num_sets());
        std::iota(set_indices.begin(), set_indices.end(), 0);
        std::sort(set_indices.begin(), set_indices.end(), [&](size_t a, size_t b)
                  { return set_cover.get_set_cost(a) < set_cover.get_set_cost(b); });
        bool solved = false;
        size_t i = 0;
        while (!solved && m_solution.size() < set_cover.get_num_sets())
        {
            auto [u, v, w] = m_links[set_indices[i]];
            add_set(m_covered_elements, u, v);
            solved = this->check_solved(m_covered_elements);
            m_solution.insert(set_indices[i]);
            i++;
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

// Partial specialization for ILP with SetCoverPseudo
// TODO: under construction
template <typename link_node_T, typename link_edge_T>
class SetCoverSolverILP<SetCoverPseudo<link_node_T, link_edge_T>>
    : public SetCoverSolver<SetCoverSolverILP<SetCoverPseudo<link_node_T, link_edge_T>>, SetCoverPseudo<link_node_T, link_edge_T>>
{
public:
    using Base = SetCoverSolver<SetCoverSolverILP<SetCoverPseudo<link_node_T, link_edge_T>>, SetCoverPseudo<link_node_T, link_edge_T>>;
    using Base::add_set;
    using Base::m_chosen_sets;
    using Base::m_covered_elements;
    using Base::m_feasible;
    using Base::m_solution;
    using Base::m_total_covered_elements;
    using Base::NUM_ELEMENTS;
    using Base::set_cover;
    using Base::SetCoverSolver;

    std::vector<std::tuple<size_t, size_t, double>> m_links;

    void solve()
    {
        double start = omp_get_wtime();
        auto csr = WeightedCRFGraph(set_cover.link_vertices,
                                    set_cover.link_edges,
                                    set_cover.link_weights);
        m_links = csr.csr_to_vec_links();
        double end = omp_get_wtime();

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
        // coverage constraints
        std::vector<GRBLinExpr> cover_expr = std::vector<GRBLinExpr>(set_cover.get_num_elements());
        for (size_t i = 0; i < set_cover.get_num_sets(); i++)
        {
            auto [u, v, w] = m_links[i];
            for (size_t k{0}; k < set_cover.n_cols; k++)
            {
                ull set_coverage = set_cover.min_cuts[u * set_cover.n_cols + k] ^
                                   set_cover.min_cuts[v * set_cover.n_cols + k];
                while (set_coverage)
                {
                    size_t bit_pos = std::countr_zero(set_coverage);
                    size_t elem = k * (8 * sizeof(ull)) + bit_pos;
                    // shouldn't be necessary
                    if (elem < set_cover.get_num_elements())
                    {
                        cover_expr[elem] += vars[i];
                    }
                    set_coverage &= set_coverage - 1; // clear the least significant bit set compiler will optimize (i checked the assembly)
                }
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
        if (status == GRB_OPTIMAL || status == GRB_SUBOPTIMAL)
        {
            for (size_t i = 0; i < set_cover.get_num_sets(); i++)
            {
                double val = vars[i].get(GRB_DoubleAttr_X);
                if (val > 0.5)
                    m_solution.insert(i);
            }
            m_feasible = true;
        }
        else
        {
            m_feasible = false;
        }
    }
};