#pragma once

#include "HeiConnect/set_cover/set_cover_bit.hpp"
#include "HeiConnect/set_cover/set_cover_csr.hpp"
#include "HeiConnect/set_cover/set_cover_cyc.hpp"
#include "HeiConnect/set_cover/set_cover_oracle.hpp"
#include "HeiConnect/set_cover/set_cover_pseudo.hpp"
#include "HeiConnect/set_cover/trimmer.hpp"
#include "HeiConnect/set_cover/solver_greedy_context.hpp"

/*
class SetCoverSolver<SetCoverBit>
{
public:
    using SetCoverType = SetCoverBit;

    void add_set(size_t set_index)
    {
        size_t k{0};
        size_t covered = 0;
#if __AVX512__
        for (; k + 7 < this->set_cover->n_cols; k += 8)
        {
            __m512i set_vec = _mm512_loadu_si512((__m512i *)(this->set_cover->set_cover.data() + set_index *
this->set_cover->n_cols + k));
            __m512i covered_vec = _mm512_loadu_si512((__m512i *)(&m_covered_elements[k]));
            __m512i new_bits = _mm512_andnot_si512(covered_vec, set_vec);
            for (int i = 0; i < 8; i++)
            {
                covered += std::popcount(((ull *)&new_bits)[i]);
            }
            __m512i new_coverage = _mm512_or_si512(set_vec, covered_vec);
            _mm512_storeu_si512((__m512i *)(&m_covered_elements[k]), new_coverage);
        }
#elif __AVX2__
        for (; k + 3 < this->set_cover->n_cols; k += 4)
        {
            __m256i set_vec = _mm256_loadu_si256((__m256i *)(this->set_cover->set_cover.data() + set_index *
this->set_cover->n_cols + k));
            __m256i covered_vec = _mm256_loadu_si256((__m256i *)(&m_covered_elements[k]));
            __m256i new_bits = _mm256_andnot_si256(covered_vec, set_vec);
            for (int i = 0; i < 4; i++)
            {
                covered += std::popcount(((ull *)&new_bits)[i]);
            }
            __m256i new_coverage = _mm256_or_si256(set_vec, covered_vec);
            _mm256_storeu_si256((__m256i *)(&m_covered_elements[k]), new_coverage);
        }
#endif
        for (; k < this->set_cover->n_cols; k++)
        {
            const ull set_word = this->set_cover->set_cover[set_index * this->set_cover->n_cols + k];
            covered += std::popcount(set_word & ~m_covered_elements[k]);
            m_covered_elements[k] |= set_word;
        }

        m_total_covered_elements += covered;
        m_solution.insert(set_index);
    }

protected:
    size_t cover_count(size_t set_index)
    {
        size_t n_covered_els = 0;
        size_t k{0};
#ifdef __AVX512VPOPCNTDQ__
        for (; k + 7 < this->set_cover->n_cols; k += 8)
        {
            __m512i vec_set = _mm512_loadu_si512((__m512i *)(this->set_cover->set_cover.data() + set_index *
this->set_cover->n_cols + k));
            __m512i vec_covered = _mm512_loadu_si512(&m_covered_elements[k]);
            __m512i vec_new_bits = _mm512_andnot_si512(vec_covered, vec_set);
            __m512i vec_popcnt = _mm512_popcnt_epi64(vec_new_bits);
            n_covered_els += _mm512_reduce_add_epi64(vec_popcnt);
        }
#elif __AVX512__
        for (; k + 7 < this->set_cover->n_cols; k += 8)
        {
            __m512i vec_set = _mm512_loadu_si512((__m512i *)(this->set_cover->set_cover.data() + set_index *
this->set_cover->n_cols + k));
            __m512i vec_covered = _mm512_loadu_si512(&m_covered_elements[k]);
            __m512i vec_new_bits = _mm512_andnot_si512(vec_covered, vec_set);
            for (int i = 0; i < 8; i++)
            {
                n_covered_els += std::popcount(((ull *)&vec_new_bits)[i]);
            }
        }
#elif __AVX2__
        for (; k + 3 < this->set_cover->n_cols; k += 4)
        {
            __m256i vec_set = _mm256_loadu_si256((__m256i *)(this->set_cover->set_cover.data() + set_index *
this->set_cover->n_cols + k));
            __m256i covered_vec = _mm256_loadu_si256((__m256i *)(&m_covered_elements[k]));
            __m256i new_bits = _mm256_andnot_si256(covered_vec, vec_set);
            for (int i = 0; i < 4; i++)
            {
                n_covered_els += std::popcount(((ull *)&new_bits)[i]);
            }
        }
#endif
        {
            const ull *set_data = this->set_cover->set_cover.data() + set_index * this->set_cover->n_cols;
            for (; k < this->set_cover->n_cols; k++)
            {
                ull covered_cuts = set_data[k] & ~m_covered_elements[k];
                n_covered_els += std::popcount(covered_cuts);
            }
        }
        return n_covered_els;
    }

    void on_greedy_solve_finished() {}

    bool check_solved(std::vector<ull> &covered_els)
    {
        size_t k{0};
#if __AVX512__
        for (; k + 7 < this->set_cover->n_cols; k += 8)
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
        for (; k + 3 < this->set_cover->n_cols; k += 4)
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
        for (; k < this->set_cover->n_cols; k++)
        {
            if (covered_els[k] != 0xFFFFFFFFFFFFFFFFULL)
            {
                return false;
            }
        }
        return true;
    }

    SolverOutput build_output() const
    {
        SolverOutput output;
        output.status = m_feasible ? SolverStatus::Feasible : SolverStatus::Unknown;
        output.chosen_sets = m_solution;
        output.objective_value = get_solution_cost();
        output.covered_elements = m_total_covered_elements;
        return output;
    }

protected:
    size_t NUM_ELEMENTS;
    const SetCoverType *set_cover = nullptr;
    SetCoverGreedySolver<SetCoverSolver> m_greedy;
    SolverConfig m_config{};
    SolverOutput m_output{};
    size_t m_total_covered_elements = 0;
    std::unordered_set<size_t> m_solution{};
    std::vector<ull> m_covered_elements;
    std::vector<bool> m_chosen_sets{};
    bool m_feasible = false;
};

template <typename Derived, class link_node_T, class link_edge_T>
class SetCoverSolver<Derived, SetCoverPseudo<link_node_T, link_edge_T>>
{
public:
    using SetCoverType = SetCoverPseudo<link_node_T, link_edge_T>;

    SetCoverSolver()
        : m_greedy(*this)
    {
    }

    ~SetCoverSolver() = default;

    template <auto AddSetFn = nullptr, auto CoverCountFn = nullptr>
    void greedy_solve()
    {
        this->m_greedy.template solve<AddSetFn, CoverCountFn>(*this->set_cover);
    }

    void solve(std::shared_ptr<const SetCoverType> sc, const SolverConfig &config)
    {
        set_cover = sc.get();
        m_config = config;
        NUM_ELEMENTS = set_cover->get_num_elements();
        m_covered_elements = std::vector<ull>(set_cover->n_cols, 0);
        m_covered_elements[set_cover->n_cols - 1] = 0xFFFFFFFFFFFFFFFFULL << (set_cover->get_num_elements() % (8 *
sizeof(ull))); WeightedCRFGraph<> g = WeightedCRFGraph<>{ {set_cover->link_vertices, set_cover->link_edges},
            set_cover->link_weights};
        m_links = g.csr_to_vec_links();
        static_cast<Derived *>(this)->algorithm_solve();
    }

    virtual void algorithm_solve() = 0;

    std::unordered_set<size_t> get_solution() const noexcept
    {
        return m_solution;
    }

    double get_solution_cost() const
    {
        return m_output.objective_value;
    }

    // Interface for greedy solver
    const SetCoverType &get_set_cover() const noexcept
    {
        return *set_cover;
    }

    size_t get_total_covered_elements() const noexcept
    {
        return m_total_covered_elements;
    }

    size_t get_solution_size() const noexcept
    {
        return m_solution.size();
    }

    bool can_remove(size_t set_index)
    {
        const auto &set_cover = *this->set_cover;
        auto [u, v, w] = m_links[set_index];
        for (size_t k{0}; k < this->set_cover->n_cols; k++)
        {
            ull set_coverage = this->set_cover->min_cuts[u * this->set_cover->n_cols + k] ^
                               this->set_cover->min_cuts[v * this->set_cover->n_cols + k];
            for (const auto &other_set_index : m_solution)
            {
                if (other_set_index == set_index)
                    continue;
                auto [ou, ov, ow] = m_links[other_set_index];
                set_coverage &= ~(this->set_cover->min_cuts[ou * this->set_cover->n_cols + k] ^
                                  this->set_cover->min_cuts[ov * this->set_cover->n_cols + k]);
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
        m_total_covered_elements = 0;
        m_covered_elements = std::vector<ull>(set_cover->n_cols, 0);
        m_covered_elements[set_cover->n_cols - 1] = 0xFFFFFFFFFFFFFFFFULL << (set_cover->get_num_elements() % (8 *
sizeof(ull)));
    }

    void trim_solution()
    {
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

    bool local_search(int k)
    {
        auto current_solution_cost = get_solution_cost();
        std::vector<size_t> solution_vec = std::vector<size_t>(m_solution.begin(), m_solution.end());
        auto subsets = calculate_subsets(k, solution_vec.size());
        auto original_solution{m_solution};
        for (size_t i{0}; i < subsets.size(); i += k)
        {
            for (size_t j{0}; j < k; j++)
            {
                size_t set_index = solution_vec[subsets[i + j]];
                m_solution.erase(set_index);
            }
            init_coverage();
            for (const auto &other_set_index : m_solution)
            {
                auto [ou, ov, ow] = m_links[other_set_index];
                add_set(m_covered_elements, ou, ov);
            }
            m_total_covered_elements = calculate_covered_elements();
            greedy_solve();
            auto new_solution_cost = get_solution_cost();
            constexpr double epsilon = 1e-9;
            if (new_solution_cost < current_solution_cost - epsilon)
            {
                std::cout << k << " local search:" << "Improved cost from "
                          << std::fixed << std::setprecision(10) << current_solution_cost << " to " << new_solution_cost
<< std::endl; return true;
            }
            else
            {
                m_solution = original_solution;
            }
        }
        return false;
    }

    size_t calculate_covered_elements()
    {
        size_t total_covered = 0;
        for (size_t k{0}; k < set_cover->n_cols; k++)
        {
            ull col = m_covered_elements[k];
            total_covered += std::popcount(col);
        }
        size_t remainder = set_cover->get_num_elements() % (8 * sizeof(ull));
        size_t tail = remainder ? (8 * sizeof(ull) - remainder) : 0;
        return total_covered - tail;
    }

    bool check_solved(std::vector<ull> &covered_els)
    {
        for (size_t k{0}; k < set_cover->n_cols; k++)
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
        auto [u, v, w] = m_links[set_index];
        (void)w;
        return cover_count(static_cast<ull>(u), static_cast<ull>(v), m_covered_elements);
    }

    void add_set(size_t set_index)
    {
        auto [u, v, w] = m_links[set_index];
        (void)w;
        const size_t covered = cover_count(static_cast<ull>(u), static_cast<ull>(v), m_covered_elements);
        add_set(m_covered_elements, static_cast<ull>(u), static_cast<ull>(v));
        m_total_covered_elements += covered;
        m_solution.insert(set_index);
    }

protected:
    void on_greedy_solve_finished() {}

private:
    size_t cover_count(ull u, ull v, const std::vector<ull> &covered_elements)
    {
        const ull *u_data = &this->set_cover->min_cuts[u * this->set_cover->n_cols];
        const ull *v_data = &this->set_cover->min_cuts[v * this->set_cover->n_cols];
        size_t n_covered_cuts = 0;
#ifdef __AVX512VPOPCNTDQ__
        for (size_t k{0}; k + 7 < this->set_cover->n_cols; k += 8)
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

    void add_set(std::vector<ull> &covered_els, ull u, ull v)
    {
        for (size_t k{0}; k < set_cover.n_cols; k++)
        {
            covered_els[k] |= set_cover.min_cuts[u * set_cover.n_cols + k] ^
                              set_cover.min_cuts[v * set_cover.n_cols + k];
        }
    }

protected:
    size_t NUM_ELEMENTS;
    const SetCoverType *set_cover = nullptr;
    SetCoverGreedySolver<SetCoverSolver> m_greedy;
    size_t m_total_covered_elements = 0;
    std::unordered_set<size_t> m_solution{};
    std::vector<bool> m_chosen_sets{};
    bool m_feasible = false;
    std::vector<ull> m_covered_elements;
    std::vector<std::tuple<size_t, size_t, double>> m_links;
};

template <typename Derived, class link_node_T, class link_edge_T, class link_weight_T>
class SetCoverSolver<Derived, SetCoverOracle<link_node_T, link_edge_T, link_weight_T>>
{
public:
    using SetCoverType = SetCoverOracle<link_node_T, link_edge_T, link_weight_T>;

    SetCoverSolver()
        : m_greedy(*this)
    {
    }

    ~SetCoverSolver() = default;

    template <auto AddSetFn = nullptr, auto CoverCountFn = nullptr>
    void greedy_solve()
    {
        this->m_greedy.template solve<AddSetFn, CoverCountFn>(*this->set_cover);
    }

    void solve(std::shared_ptr<const SetCoverType> sc, const SolverConfig &config)
    {
        set_cover = sc.get();
        m_config = config;
        NUM_ELEMENTS = set_cover->get_num_elements();
        m_coverage_count = std::vector<ull>(set_cover->get_num_elements(), 0);
        WeightedCRFGraph<> g = WeightedCRFGraph<>{
            {set_cover->link_vertices, set_cover->link_edges},
            set_cover->link_weights};
        m_links = g.csr_to_vec_links();
        static_cast<Derived *>(this)->algorithm_solve();
    }

    // Virtual algorithm-specific solve (to be implemented by derived class)
    virtual void algorithm_solve() = 0;

    std::unordered_set<size_t> get_solution() const noexcept
    {
        return m_solution;
    }

    double get_solution_cost() const
    {
        return m_output.objective_value;
    }

    // Interface for greedy solver
    const SetCoverType &get_set_cover() const noexcept
    {
        return *set_cover;
    }

    size_t get_total_covered_elements() const noexcept
    {
        return m_total_covered_elements;
    }

    size_t get_solution_size() const noexcept
    {
        return m_solution.size();
    }

protected:
    void on_greedy_solve_finished() {}

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

    void trim_solution()
    {
        std::vector<size_t> solution_vec(m_solution.begin(), m_solution.end());
        std::sort(solution_vec.begin(), solution_vec.end(), [&](size_t a, size_t b)
                  { return set_cover.get_set_cost(a) > set_cover.get_set_cost(b); });
        for (const auto &set_index : solution_vec)
        {
            if (this->can_remove(set_index))
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
        const auto &covered_elements = get_or_build_set_coverage(set_index);
        for (const auto cut_idx : covered_elements)
        {
            if (m_coverage_count[cut_idx] == 0)
            {
                covered++;
            }
        }
        return covered;
    }

    void add_set(size_t set_index)
    {
        const auto &covered_elements = get_or_build_set_coverage(set_index);
        for (const auto cut_idx : covered_elements)
        {
            const size_t cov = m_coverage_count[cut_idx]++;
            if (cov == 0)
            {
                m_total_covered_elements++;
            }
        }
        m_solution.insert(set_index);
    }

    void delete_set(size_t set_index)
    {
        const auto &covered_elements = get_or_build_set_coverage(set_index);
        for (const auto cut_idx : covered_elements)
        {
            const size_t cov = --m_coverage_count[cut_idx];
            if (cov == 0)
            {
                m_total_covered_elements--;
            }
        }
        m_solution.erase(set_index);
    }

private:
    const std::vector<size_t> &get_or_build_set_coverage(size_t set_index)
    {
        auto cached = m_set_coverage_cache.find(set_index);
        if (cached != m_set_coverage_cache.end())
        {
            return cached->second;
        }

        std::vector<size_t> covered_elements;
        covered_elements.reserve(set_cover.get_num_tree_cuts() + set_cover.get_num_cycle_cuts());
        for (size_t k{0}; k < set_cover.get_num_tree_cuts(); k++)
        {
            if (set_cover.covers_tree(set_index, k))
            {
                covered_elements.emplace_back(k);
            }
        }
        for (size_t k{0}; k < set_cover.get_num_cycle_cuts(); k++)
        {
            if (set_cover.covers_cycle(set_index, k))
            {
                covered_elements.emplace_back(set_cover.get_num_tree_cuts() + k);
            }
        }
        auto [it, _] = m_set_coverage_cache.emplace(set_index, std::move(covered_elements));
        return it->second;
    }

protected:
    void on_greedy_solve_finished() {}

protected:
    size_t NUM_ELEMENTS;
    const SetCoverType *set_cover = nullptr;
    SetCoverGreedySolver<SetCoverSolver> m_greedy;
    SolverConfig m_config{};
    SolverOutput m_output{};
    size_t m_total_covered_elements = 0;
    std::unordered_set<size_t> m_solution{};
    std::vector<bool> m_chosen_sets{};
    bool m_feasible = false;
    std::vector<ull> m_covered_elements;
    std::vector<ull> m_coverage_count;
    std::vector<std::tuple<size_t, size_t, double>> m_links;
    std::unordered_map<size_t, std::vector<size_t>> m_set_coverage_cache;
};

template <typename Derived, class link_node_T, class link_edge_T, class link_weight_T>
class SetCoverSolver<Derived, SetCoverCyc<link_node_T, link_edge_T, link_weight_T>>
{
public:
    using SetCoverType = SetCoverCyc<link_node_T, link_edge_T, link_weight_T>;
    struct ArcEquivClass
    {
        cycle_pos_T cls;
        cycle_pos_T start;
        cycle_pos_T end;
    };

    SetCoverSolver()
        : m_greedy(*this)
    {
    }

    ~SetCoverSolver() = default;

    template <auto AddSetFn = nullptr, auto CoverCountFn = nullptr>
    void greedy_solve()
    {
        this->m_greedy.template solve<AddSetFn, CoverCountFn>(*this->set_cover);
    }

    void solve(std::shared_ptr<const SetCoverType> sc, const SolverConfig &config)
    {
        set_cover = sc.get();
        m_config = config;
        NUM_ELEMENTS = set_cover.get_num_elements();
        m_coverage_count = std::vector<ull>(set_cover.get_num_elements(), 0);
        WeightedCRFGraph<> g = WeightedCRFGraph<>{
            {set_cover.link_vertices, set_cover.link_edges},
            set_cover.link_weights};
        m_links = g.csr_to_vec_links();
        m_arc_equiv_classes = std::vector<std::vector<ArcEquivClass>>(set_cover.cycle_sizes.size());
        for (size_t c{0}; c < set_cover.cycle_sizes.size(); c++)
        {
            m_arc_equiv_classes[c].push_back({0, 0, static_cast<cycle_pos_T>(set_cover.cycle_sizes[c])});
            m_class_sizes.emplace_back(std::vector<cycle_pos_T>{static_cast<cycle_pos_T>(set_cover.cycle_sizes[c])});
        }
        m_class_intersects = std::vector<size_t>(*std::max_element(set_cover.cycle_sizes.begin(),
set_cover.cycle_sizes.end()), 0); static_cast<Derived *>(this)->algorithm_solve();
    }

    virtual void algorithm_solve() = 0;

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

    // Interface for greedy solver
    const SetCoverType &get_set_cover() const noexcept
    {
        return *set_cover;
    }

    size_t get_total_covered_elements() const noexcept
    {
        return m_total_covered_elements;
    }

    size_t get_solution_size() const noexcept
    {
        return m_solution.size();
    }

protected:
    void init_coverage()
    {
        m_total_covered_elements = 0;
        m_covered_elements = std::vector<ull>(set_cover.n_cols, 0);
        m_covered_elements[set_cover.n_cols - 1] = 0xFFFFFFFFFFFFFFFFULL << (set_cover.get_num_tree_cuts() % (8 *
sizeof(ull)));
    }

    void trim_solution()
    {
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

    void on_greedy_solve_finished()
    {
        const double denom = static_cast<double>(m_solution.size()) * static_cast<double>(set_cover.get_num_sets());
        const double spared_coverage = denom > 0.0 ? (1.0 -
(static_cast<double>(m_cover_count_metrics.cover_count_counter) / denom)) : 0.0; std::cout << "Cover count spare ratio:
" << spared_coverage << std::endl;
    }

private:
    struct CoverCountMetrics
    {
        size_t cover_count_counter = 0;
        double time_cover_count_total = 0.0;
        double time_cover_count_cycle_part = 0.0;
        double time_cover_count_tree_part = 0.0;
        double time_class_intersects = 0.0;
        double time_alloc_class_intersects = 0.0;
        double time_cover_count_class_accumulation = 0.0;

        double time_add_set_total = 0.0;
        double time_add_set_cycle_refine = 0.0;
        double time_add_set_insert_solution = 0.0;
        double time_add_set_update_total_covered = 0.0;
    } m_cover_count_metrics;

    std::vector<size_t> m_class_intersects;
    SolverConfig m_config{};
    SolverOutput m_output{};
    size_t m_total_covered_elements = 0;
    std::unordered_set<size_t> m_solution{};
    std::vector<bool> m_chosen_sets{};
    bool m_feasible = false;
    std::vector<ull> m_covered_elements;
    std::vector<ull> m_coverage_count;
    std::vector<std::tuple<size_t, size_t, double>> m_links;
    std::unordered_map<size_t, std::vector<size_t>> m_set_coverage_cache;

public:
    size_t cover_count_1(size_t set_index)
    {
        const double total_start = omp_get_wtime();
        m_cover_count_metrics.cover_count_counter++;

        auto [u, v, w] = m_links[set_index];
        size_t covered = 0;

        const double cycle_part_start = omp_get_wtime();
        const auto &ccs = set_cover.m_cycle_crosses[set_index];
        for (const auto &cc : ccs)
        {
            auto [cycle, a, b] = cc;
            std::vector<ArcEquivClass> &arcs = m_arc_equiv_classes[cycle];
            double section_start = omp_get_wtime();
            m_class_intersects.assign(m_class_sizes[cycle].size(), 0);
            m_cover_count_metrics.time_alloc_class_intersects += omp_get_wtime() - section_start;
            section_start = omp_get_wtime();
            for (size_t arc_i{0}; arc_i < arcs.size(); arc_i++)
            m_runtime_seconds = end - start;
            m_objective_value = get_solution_cost();
            m_covered_elements = count_covered_elements();
            m_status = (m_covered_elements == m_set_cover->get_num_elements()) ? SolverStatus::Feasible :
SolverStatus::Unknown;
                {
                    m_class_intersects[cls] += std::max(cycle_pos_T{0}, intersect_end - intersect_start);
        const SolutionType &get_solution() const
                else
            return m_solution;
                    if (intersect_start < intersect_end)
                    {
        double get_solution_cost() const
                    }
            if (!m_set_cover)
            {
                return 0.0;
            }

            double total_cost = 0.0;
            for (const auto set_index : m_solution.get_solution())
            {
                total_cost += m_set_cover->get_set_cost(set_index);
            }
            return total_cost;
            }
            m_cover_count_metrics.time_class_intersects += omp_get_wtime() - section_start;
        size_t get_covered_elements() const noexcept
            section_start = omp_get_wtime();
            return m_covered_elements;
        }

        double get_runtime_seconds() const noexcept
        {
            return m_runtime_seconds;
        }

        SolverStatus get_status() const noexcept
        {
            return m_status;
            ull coverage = set_cover.m_tree_partition_matrix[u * set_cover.get_num_tree_cuts() + k] ^
set_cover.m_tree_partition_matrix[v * set_cover.get_num_tree_cuts() + k]; covered += std::popcount(coverage &
~m_covered_elements[k]);
        }
        size_t count_covered_elements() const
        {
            if (!m_set_cover)
            {
                return 0;
            }

            std::vector<bool> covered(m_set_cover->get_num_elements(), false);
            for (const auto set_index : m_solution.get_solution())
            {
                m_set_cover->forEachElement(set_index, [&](size_t element)
                {
                    covered[element] = true;
                });
            }
            return std::count(covered.begin(), covered.end(), true);
        }

        m_cover_count_metrics.time_cover_count_tree_part += omp_get_wtime() - tree_part_start;
        return covered;
    }
        size_t m_covered_elements = 0;
        double m_runtime_seconds = 0.0;
        double m_objective_value = 0.0;
        SolverStatus m_status = SolverStatus::Unknown;

private:
    static constexpr cycle_pos_T INVALID_CLASS = std::numeric_limits<cycle_pos_T>::max();

    static void merge_adjacent_arcs(std::vector<ArcEquivClass> &arcs)
    {
        if (arcs.empty())
            return;

        std::vector<ArcEquivClass> merged;
        merged.reserve(arcs.size());
        merged.push_back(arcs[0]);

        for (size_t i = 1; i < arcs.size(); ++i)
        {
            auto &back = merged.back();
            if (back.cls == arcs[i].cls && back.end == arcs[i].start)
            {
                back.end = arcs[i].end;
            }
            else
            {
                merged.push_back(arcs[i]);
            }
        }
        arcs.swap(merged);
    }

    size_t refine_cycle_partition(size_t cycle, cycle_pos_T a, cycle_pos_T b)
    {
        auto &arcs = m_arc_equiv_classes[cycle];
        auto &class_sizes = m_class_sizes[cycle];

        const size_t old_num_classes = class_sizes.size();

        std::vector<cycle_pos_T> inside(old_num_classes, 0);

        for (const auto &arc : arcs)
        {
            const cycle_pos_T l = std::max(arc.start, a);
            const cycle_pos_T r = std::min(arc.end, b);
            if constexpr (std::signed_integral<cycle_pos_T>)
            {
                inside[arc.cls] += std::max(cycle_pos_T{0}, r - l);
            }
            else
            {
                if (l < r)
                {
                    inside[arc.cls] += (r - l);
                }
            }
        }

        std::vector<cycle_pos_T> split_to(old_num_classes, INVALID_CLASS);

        size_t gained = 0;
        for (size_t c = 0; c < old_num_classes; ++c)
        {
            const cycle_pos_T in = inside[c];
            const cycle_pos_T out = class_sizes[c] - in;

            gained += static_cast<size_t>(in) * static_cast<size_t>(out);

            if (in > 0 && out > 0)
            {
                const cycle_pos_T new_cls = static_cast<cycle_pos_T>(class_sizes.size());

                split_to[c] = new_cls;

                class_sizes[c] = out;
                class_sizes.push_back(in);
            }
        }

        std::vector<ArcEquivClass> new_arcs;
        new_arcs.reserve(arcs.size() * 2 + 4);

        for (const auto &arc : arcs)
        {
            const cycle_pos_T l = std::max(arc.start, a);
            const cycle_pos_T r = std::min(arc.end, b);

            if (!(l < r))
            {
                new_arcs.push_back(arc);
                continue;
            }

            const cycle_pos_T new_cls = split_to[arc.cls];

            if (new_cls == INVALID_CLASS)
            {
                new_arcs.push_back(arc);
                continue;
            }

            if (arc.start < l)
            {
                new_arcs.push_back({arc.cls, arc.start, l});
            }

            new_arcs.push_back({new_cls, l, r});

            if (r < arc.end)
            {
                new_arcs.push_back({arc.cls, r, arc.end});
            }
        }

        arcs.swap(new_arcs);
        merge_adjacent_arcs(arcs);

        return gained;
    }

public:
    void add_set_1(size_t set_index)
    {
        const double total_start = omp_get_wtime();

        auto [u, v, w] = m_links[set_index];

        size_t gained = 0;

        const double cycle_refine_start = omp_get_wtime();
        const auto &ccs = set_cover.m_cycle_crosses[set_index];
        for (const auto &cc : ccs)
        {
            auto [cycle, a, b] = cc;
            gained += refine_cycle_partition(cycle, a, b);
        }
        m_cover_count_metrics.time_add_set_cycle_refine += omp_get_wtime() - cycle_refine_start;

        double section_start = omp_get_wtime();
        m_solution.insert(set_index);
        m_cover_count_metrics.time_add_set_insert_solution += omp_get_wtime() - section_start;

        section_start = omp_get_wtime();
        m_total_covered_elements += gained;
        m_cover_count_metrics.time_add_set_update_total_covered += omp_get_wtime() - section_start;

        m_cover_count_metrics.time_add_set_total += omp_get_wtime() - total_start;
    }

    size_t cover_count(size_t set_index)
    {
        m_cover_count_metrics.cover_count_counter++;
        size_t covered = 0;
        auto [u, v, w] = m_links[set_index];
        for (size_t k{0}; k < set_cover.get_num_tree_cuts(); k++)
        {
            ull coverage = set_cover.m_tree_partition_matrix[u * set_cover.n_cols + k] ^
                           set_cover.m_tree_partition_matrix[v * set_cover.n_cols + k];
            covered += std::popcount(coverage & ~m_covered_elements[k]);
        }
        const auto &cc = set_cover.m_cycle_crosses[set_index];
        for (size_t i{0}; i < cc.size(); i++)
        {
            auto [c, a, b] = cc[i];
            auto c_size = set_cover.cycle_sizes[c];
            const size_t a_pos = static_cast<size_t>(a);
            const size_t b_pos = static_cast<size_t>(b);
            for (size_t j{0}; j < a_pos; j++)
            {
                for (size_t k = a_pos; k < b_pos; k++)
                {
                    if (set_cover.m_cycle_coverages[c][j * c_size + k] == 0)
                    {
                        covered++;
                    }
                }
            }
            for (size_t j = a_pos; j < b_pos; j++)
            {
                for (size_t k = b_pos; k < c_size; k++)
                {
                    if (set_cover.m_cycle_coverages[c][j * c_size + k] == 0)
                    {
                        covered++;
                    }
                }
            }
        }

        return covered;
    }

    void add_set(size_t set_index)
    {
        auto [u, v, w] = m_links[set_index];
        auto covered = 0ull;
        for (size_t k{0}; k < set_cover.n_cols; k++)
        {
            auto cov = set_cover.m_tree_partition_matrix[u * set_cover.n_cols + k] ^ set_cover.m_tree_partition_matrix[v
* set_cover.n_cols + k]; covered += std::popcount(cov & ~m_covered_elements[k]); m_covered_elements[k] |= cov;
        }
        const auto &cc = set_cover.m_cycle_crosses[set_index];
        for (size_t i{0}; i < cc.size(); i++)
        {
            auto [c, a, b] = cc[i];
            auto c_size = set_cover.cycle_sizes[c];
            const size_t a_pos = static_cast<size_t>(a);
            const size_t b_pos = static_cast<size_t>(b);
            for (size_t j{0}; j < a_pos; j++)
            {
                for (size_t k = a_pos; k < b_pos; k++)
                {
                    if (set_cover.m_cycle_coverages[c][j * c_size + k]++ == 0)
                    {
                        covered++;
                    }
                }
            }
            for (size_t j = a_pos; j < b_pos; j++)
            {
                for (size_t k = b_pos; k < c_size; k++)
                {
                    if (set_cover.m_cycle_coverages[c][j * c_size + k]++ == 0)
                    {
                        covered++;
                    }
                }
            }
        }
        m_total_covered_elements += covered;
        m_solution.insert(set_index);
    }

    void delete_set(size_t set_index)
    {
        (void)set_index;
        return;
    }

    bool can_remove(size_t set_index)
    {
        (void)set_index;
        return true;
    }

    void print_metrics() const
    {
        const double cover_total = m_cover_count_metrics.time_cover_count_total;
        const double add_total = m_cover_count_metrics.time_add_set_total;
        auto pct = [](double part, double total) -> double
        {
            return total > 0.0 ? (100.0 * part / total) : 0.0;
        };

        std::cout << "Cover count calls: " << m_cover_count_metrics.cover_count_counter << std::endl;

        std::cout << "Cover count total time: " << cover_total << " seconds" << std::endl;
        std::cout << "  Cover count cycle part time: " << m_cover_count_metrics.time_cover_count_cycle_part << " seconds
(" << pct(m_cover_count_metrics.time_cover_count_cycle_part, cover_total) << "%)" << std::endl; std::cout << "  Cover
count tree part time: " << m_cover_count_metrics.time_cover_count_tree_part << " seconds (" <<
pct(m_cover_count_metrics.time_cover_count_tree_part, cover_total) << "%)" << std::endl; std::cout << "  Cover count
class intersects time: " << m_cover_count_metrics.time_class_intersects << " seconds (" <<
pct(m_cover_count_metrics.time_class_intersects, cover_total) << "%)" << std::endl; std::cout << "  Cover count allocate
class intersects time: " << m_cover_count_metrics.time_alloc_class_intersects << " seconds (" <<
pct(m_cover_count_metrics.time_alloc_class_intersects, cover_total) << "%)" << std::endl; std::cout << "  Cover count
class accumulation time: " << m_cover_count_metrics.time_cover_count_class_accumulation << " seconds (" <<
pct(m_cover_count_metrics.time_cover_count_class_accumulation, cover_total) << "%)" << std::endl;

        std::cout << "Add set total time: " << add_total << " seconds" << std::endl;
        std::cout << "  Add set cycle refine time: " << m_cover_count_metrics.time_add_set_cycle_refine << " seconds ("
<< pct(m_cover_count_metrics.time_add_set_cycle_refine, add_total) << "%)" << std::endl; std::cout << "  Add set insert
solution time: " << m_cover_count_metrics.time_add_set_insert_solution << " seconds (" <<
pct(m_cover_count_metrics.time_add_set_insert_solution, add_total) << "%)" << std::endl; std::cout << "  Add set
covered-elements update time: " << m_cover_count_metrics.time_add_set_update_total_covered << " seconds (" <<
pct(m_cover_count_metrics.time_add_set_update_total_covered, add_total) << "%)" << std::endl;
    }

protected:
    size_t NUM_ELEMENTS;
    const SetCoverType *set_cover = nullptr;
    SetCoverGreedySolver<SetCoverSolver> m_greedy;
    size_t m_total_covered_elements = 0;
    std::unordered_set<size_t> m_solution{};
    std::vector<bool> m_chosen_sets{};
    bool m_feasible = false;
    std::vector<ull> m_covered_elements;
    std::vector<ull> m_coverage_count;
    std::vector<std::tuple<size_t, size_t, double>> m_links;
    std::vector<std::vector<ArcEquivClass>> m_arc_equiv_classes;
    std::vector<std::vector<cycle_pos_T>> m_class_sizes;
};
*/
