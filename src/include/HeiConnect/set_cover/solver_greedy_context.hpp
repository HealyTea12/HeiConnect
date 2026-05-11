#pragma once

#include "HeiConnect/set_cover/common.hpp"

#include <memory>
#include <vector>

template <typename SetCoverType, typename Solution>
    requires SetCoverCon<SetCoverType> && ForEachElementCon<SetCoverType>
class BasicContext
{
    using CoverCount = int;

public:
    explicit BasicContext(std::shared_ptr<const SetCoverType> set_cover)
        : m_setCover(std::move(set_cover)), m_coveredCount(m_setCover->get_num_elements(), 0)
    {
    }

    void add_set(size_t set_index)
    {
        m_setCover->forEachElement(set_index, [this](typename SetCoverType::ElementID element)
                                   { ++m_coveredCount[element]; });
    }

    void remove_set(size_t set_index)
    {
        m_setCover->forEachElement(set_index, [this](typename SetCoverType::ElementID element)
                                   { --m_coveredCount[element]; });
    }

    size_t cover_count(size_t set_index) const
    {
        size_t count = 0;
        m_setCover->forEachElement(set_index, [this, &count](typename SetCoverType::ElementID element)
                                   {
            if (m_coveredCount[element] == 0)
            {
                ++count;
            } });
        return count;
    }

    bool can_remove(size_t set_index, const Solution &) const
    {
        bool removable = true;
        m_setCover->forEachElement(set_index, [this, &removable](typename SetCoverType::ElementID element)
                                   {
            if (m_coveredCount[element] == 1)
            {
                removable = false;
            } });
        return removable;
    }

    bool is_element_covered(size_t element_index) const
    {
        return m_coveredCount[element_index] > 0;
    }

    size_t get_total_covered_elements() const
    {
        return std::count_if(m_coveredCount.begin(), m_coveredCount.end(), [](CoverCount count)
                             { return count > 0; });
    }

    void reset()
    {
        std::fill(m_coveredCount.begin(), m_coveredCount.end(), 0);
    }

private:
    std::shared_ptr<const SetCoverType> m_setCover;
    std::vector<CoverCount> m_coveredCount;
};

// Specialized context for bit-packed representation (SetCoverBit)
template <typename SetCoverType, typename Solution>
class BitPackedContext
{
    using CoverCount = int;
    using ull = unsigned long long; // TODO: should actually be the word size of the underlying type
    constexpr static size_t WORD_BITS = 8 * sizeof(ull);

public:
    explicit BitPackedContext(std::shared_ptr<const SetCoverType> set_cover)
        : m_setCover(std::move(set_cover)), m_coverageMask(m_setCover->get_n_cols(), 0)
    {
    }

    void add_set(typename SetCoverType::SetID set_index)
    {
        size_t col = 0;
        m_setCover->forEachElementBitMasked(set_index, [this, &col](size_t word)
                                            {
            m_coverageMask[col] |= word;
            col++; });
    }

    // ATTENTION: Need to be very careful when calling this one because it breaks the underlying state
    // so subsequent calls to add set and cover count will be incorrect
    // this method is required so that we can pass the same context to the trimmer
    void remove_set(typename SetCoverType::SetID set_index)
    {
        // Does nothing
        return;
    }

    size_t cover_count(SetCoverType::SetID set_index) const
    {
        size_t count = 0;
        size_t col = 0;
        m_setCover->forEachElementBitMasked(set_index, [&](size_t word)
                                            {
            ull new_bits = word & ~m_coverageMask[col];
            count += std::popcount(new_bits);
            col++; });
        return count;
    }

    // this repeats the same as the trimmer context for bitpacked maybe remove repetition
    bool can_remove(SetCoverType::SetID set_index, const Solution &solution) const
    {
        const auto &solSet = solution.get_solution();
        for (size_t k{0}; k < m_setCover->get_n_cols(); k++)
        {
            ull set_coverage = m_setCover->get_col(set_index, k);

            for (const auto &other_set_index : solSet)
            {
                if (other_set_index == set_index)
                    continue;
                set_coverage &= ~(m_setCover->get_col(other_set_index, k));
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

    size_t get_total_covered_elements() const
    {
        size_t total = 0;
        for (size_t col = 0; col < m_setCover->get_n_cols(); col++)
        {
            total += std::popcount(m_coverageMask[col]);
        }
        return total;
    }

    bool is_element_covered(size_t element_index) const
    {
        size_t word_index = element_index / WORD_BITS;
        size_t bit_index = element_index % WORD_BITS;
        return (m_coverageMask[word_index] & (1ULL << bit_index)) != 0;
    }

    void reset()
    {
        std::fill(m_coverageMask.begin(), m_coverageMask.end(), 0);
    }

private:
    std::shared_ptr<const SetCoverType> m_setCover;
    std::vector<ull> m_coverageMask;
};

template <typename SetCoverType1, typename SetCoverType2,
          typename FirstContext, typename SecondContext, typename SolutionType>
class ContextDouble
{
public:
    explicit ContextDouble(std::shared_ptr<const SetCoverDouble<SetCoverType1, SetCoverType2>> set_cover)
        : m_context1(std::make_shared<FirstContext>(std::make_shared<SetCoverType1>(set_cover->first()))),
          m_context2(std::make_shared<SecondContext>(std::make_shared<SetCoverType2>(set_cover->second()))),
          m_setCover(set_cover)
    {
    }

    void add_set(size_t set_index)
    {
        m_context1->add_set(set_index);
        m_context2->add_set(set_index);
    }

    void remove_set(size_t set_index)
    {
        m_context1->remove_set(set_index);
        m_context2->remove_set(set_index);
    }

    bool can_remove(size_t set_index, const SolutionType &solution) const
    {
        return m_context1->can_remove(set_index, solution) && m_context2->can_remove(set_index, solution);
    }

    size_t cover_count(size_t set_index) const
    {
        return m_context1->cover_count(set_index) + m_context2->cover_count(set_index);
    }

    size_t get_total_covered_elements() const
    {
        return m_context1->get_total_covered_elements() + m_context2->get_total_covered_elements();
    }

    // this might be incorrect
    bool is_element_covered(size_t element_index) const
    {
        size_t firstElements = m_setCover->first().get_num_elements();
        if (element_index < firstElements)
        {
            return m_context1->is_element_covered(element_index);
        }
        else
        {
            return m_context2->is_element_covered(element_index - firstElements);
        }
    }

    void reset()
    {
        m_context1->reset();
        m_context2->reset();
    }

private:
    std::shared_ptr<FirstContext> m_context1;
    std::shared_ptr<SecondContext> m_context2;
    std::shared_ptr<const SetCoverDouble<SetCoverType1, SetCoverType2>> m_setCover;
};

// specialization for set cover cyc
template <class link_node_T, class link_edge_T, class link_weight_T, typename SolutionType>
class CycContext
{
    using CycPos = size_t;
    struct ArcEquivClass
    {
        CycPos cls;
        CycPos start;
        CycPos end;
    };

public:
    CycContext(std::shared_ptr<SetCoverCyc<link_node_T, link_edge_T, link_weight_T>> sc)
        : m_setCover(std::move(sc))
    {
        const auto &cycleSizes = m_setCover->get_cycle_sizes();
        m_arcEquivClasses = std::vector<std::vector<ArcEquivClass>>(
            cycleSizes.size());
        for (size_t c{0}; c < cycleSizes.size(); c++)
        {
            m_arcEquivClasses[c].push_back(ArcEquivClass{0, 0, static_cast<CycPos>(cycleSizes[c])});
            m_classSizes.emplace_back(std::vector<CycPos>{static_cast<CycPos>(cycleSizes[c])});
        }
        m_classIntersects = std::vector<size_t>(*std::max_element(cycleSizes.begin(), cycleSizes.end()), 0);
    }

    size_t cover_count(size_t set_index)
    {
        size_t covered = 0;
        const auto &ccs = m_setCover->get_cycle_crosses()[set_index];
        for (const auto &cc : ccs)
        {
            auto [cycle, a, b] = cc;
            std::vector<ArcEquivClass> &arcs = m_arcEquivClasses[cycle];
            m_classIntersects.assign(m_classSizes[cycle].size(), 0);
            for (size_t arc_i{0}; arc_i < arcs.size(); arc_i++)
            {
                auto [cls, start, end] = arcs[arc_i];
                auto intersect_start = std::max(start, static_cast<CycPos>(a));
                auto intersect_end = std::min(end, static_cast<CycPos>(b));
                if constexpr (std::signed_integral<CycPos>)
                {
                    m_classIntersects[cls] += std::max(CycPos{0}, intersect_end - intersect_start);
                }
                else
                {
                    if (intersect_start < intersect_end)
                    {
                        m_classIntersects[cls] += intersect_end - intersect_start;
                    }
                }
            }

            for (size_t i{0}; i < m_classSizes[cycle].size(); i++)
            {
                covered += m_classIntersects[i] * (m_classSizes[cycle][i] - m_classIntersects[i]);
            }
        }

        return covered;
    }

    void add_set(size_t set_index)
    {
        size_t gained = 0;

        const auto &ccs = m_setCover->get_cycle_crosses()[set_index];
        for (const auto &cc : ccs)
        {
            auto [cycle, a, b] = cc;
            gained += refine_cycle_partition(cycle, a, b);
        }

        m_totalCoveredElements += gained;
    }

    // TODO: right now there is no way of going backwards
    // or calculating if a set is removable in this context
    // bool can_remove(size_t) const
    //{
    //    return false;
    //}

    // void remove_set(size_t)
    //{
    // }

    size_t get_total_covered_elements() const
    {
        return m_totalCoveredElements;
    }

private:
    static constexpr CycPos INVALID_CLASS = std::numeric_limits<CycPos>::max();

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

    size_t refine_cycle_partition(size_t cycle, CycPos a, CycPos b)
    {
        auto &arcs = m_arcEquivClasses[cycle];
        auto &class_sizes = m_classSizes[cycle];

        const size_t old_num_classes = class_sizes.size();

        std::vector<CycPos> inside(old_num_classes, 0);

        for (const auto &arc : arcs)
        {
            const CycPos l = std::max(arc.start, a);
            const CycPos r = std::min(arc.end, b);
            if constexpr (std::signed_integral<CycPos>)
            {
                inside[arc.cls] += std::max(CycPos{0}, r - l);
            }
            else
            {
                if (l < r)
                {
                    inside[arc.cls] += (r - l);
                }
            }
        }

        std::vector<CycPos> split_to(old_num_classes, INVALID_CLASS);

        size_t gained = 0;
        for (size_t c = 0; c < old_num_classes; ++c)
        {
            const CycPos in = inside[c];
            const CycPos out = class_sizes[c] - in;

            gained += static_cast<size_t>(in) * static_cast<size_t>(out);

            if (in > 0 && out > 0)
            {
                const CycPos new_cls = static_cast<CycPos>(class_sizes.size());

                split_to[c] = new_cls;

                class_sizes[c] = out;
                class_sizes.emplace_back(in);
            }
        }

        std::vector<ArcEquivClass> new_arcs;
        new_arcs.reserve(arcs.size() * 2 + 4);

        for (const auto &arc : arcs)
        {
            const CycPos l = std::max(arc.start, a);
            const CycPos r = std::min(arc.end, b);

            if (!(l < r))
            {
                new_arcs.push_back(arc);
                continue;
            }

            const CycPos new_cls = split_to[arc.cls];

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

private:
    std::shared_ptr<SetCoverCyc<link_node_T, link_edge_T, link_weight_T>> m_setCover;

    // this vector is used to store intersections during cover count cacluatlions to avoid multiple allocs
    std::vector<size_t> m_classIntersects;
    std::vector<std::vector<ArcEquivClass>> m_arcEquivClasses;
    std::vector<std::vector<CycPos>> m_classSizes;
    size_t m_totalCoveredElements = 0;
};