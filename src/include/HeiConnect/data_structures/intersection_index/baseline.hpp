#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

using IntersectionInterval = std::tuple<size_t, size_t>;

inline std::ptrdiff_t interval_overlap(IntersectionInterval left, IntersectionInterval right)
{
    const auto [a, b] = left;
    const auto [c, d] = right;
    const size_t overlap_begin = std::max(a, c);
    const size_t overlap_end = std::min(b, d);
    if (overlap_end < overlap_begin)
    {
        return -1;
    }
    return static_cast<std::ptrdiff_t>(overlap_end - overlap_begin);
}

inline bool intervals_touch_or_cross(IntersectionInterval left, IntersectionInterval right)
{
    const auto [a, b] = left;
    const auto [c, d] = right;
    const bool share_endpoint = a == c || a == d || b == c || b == d;
    const std::ptrdiff_t overlap = interval_overlap(left, right);
    const size_t smaller_interval = std::min(b - a, d - c);
    const bool cross = overlap > 0 && static_cast<size_t>(overlap) < smaller_interval;
    return share_endpoint || cross;
}

struct IntersectionIndexMetrics
{
    size_t queries{};
    size_t callbacks{};
    size_t intervals_popped{};
    size_t candidates_inspected{};
    size_t subtrees_pruned_by_level{};

    void add(const IntersectionIndexMetrics& other)
    {
        queries += other.queries;
        callbacks += other.callbacks;
        intervals_popped += other.intervals_popped;
        candidates_inspected += other.candidates_inspected;
        subtrees_pruned_by_level += other.subtrees_pruned_by_level;
    }
};

struct IntersectionRecord
{
    IntersectionInterval interval;
    size_t level;
};

template<int RecordStatsLevel = 0>
class BaseIntersectionIdx
{
public:
    virtual ~BaseIntersectionIdx() = default;
    using Interval = IntersectionInterval;
    using Index = size_t;

    virtual std::unique_ptr<BaseIntersectionIdx<RecordStatsLevel>> make(
        const std::vector<IntersectionRecord>& records) const = 0;

    virtual std::unique_ptr<BaseIntersectionIdx<RecordStatsLevel>> makeEmpty(
        const std::vector<IntersectionRecord>& possible_records) const = 0;

    virtual Index addInterval(IntersectionRecord record) = 0;

    virtual void forEachIntersection(
        std::function<void(Index, Interval)> callback,
        Interval query,
        size_t exclusive_level) const = 0;

    virtual void popInterval(Interval interval) = 0;
    virtual void popInterval(Index idx) = 0;
    virtual void activateInterval(Index idx, size_t level) = 0;
    virtual void clear() = 0;
    virtual void reset() = 0;
    virtual IntersectionIndexMetrics emit_metrics() const = 0;
};

template<int RecordStatsLevel = 0>
class BaselineIntersectionIdx : public BaseIntersectionIdx<RecordStatsLevel>
{
public:
    using typename BaseIntersectionIdx<RecordStatsLevel>::Index;
    using typename BaseIntersectionIdx<RecordStatsLevel>::Interval;

    BaselineIntersectionIdx() = default;

    explicit BaselineIntersectionIdx(const std::vector<IntersectionRecord>& records) : m_records(records)
    {
        m_active.resize(records.size(), true);
        m_initial_levels.reserve(records.size());
        for (const IntersectionRecord& record : records)
        {
            m_initial_levels.push_back(record.level);
        }
    }

    explicit BaselineIntersectionIdx(size_t capacity) : m_incremental(true)
    {
        m_records.reserve(capacity);
        m_active.reserve(capacity);
        m_initial_levels.reserve(capacity);
    }

    std::unique_ptr<BaseIntersectionIdx<RecordStatsLevel>> make(
        const std::vector<IntersectionRecord>& records) const override
    {
        return std::make_unique<BaselineIntersectionIdx<RecordStatsLevel>>(records);
    }

    std::unique_ptr<BaseIntersectionIdx<RecordStatsLevel>> makeEmpty(
        const std::vector<IntersectionRecord>& possible_records) const override
    {
        return std::make_unique<BaselineIntersectionIdx<RecordStatsLevel>>(possible_records.size());
    }

    Index addInterval(IntersectionRecord record) override
    {
        if (!m_records.empty() && m_records.back().level > record.level)
        {
            throw std::invalid_argument("Intervals must be added in nondecreasing level order.");
        }

        const Index index = m_records.size();
        m_records.push_back(record);
        m_active.push_back(true);
        m_initial_levels.push_back(record.level);
        return index;
    }

    void forEachIntersection(
        std::function<void(Index, Interval)> callback,
        Interval query,
        size_t exclusive_level) const override
    {
        if constexpr (RecordStatsLevel > 0)
        {
            m_metrics.queries++;
        }
        for (size_t i{0}; i < m_records.size(); ++i)
        {
            if (m_incremental && m_records[i].level >= exclusive_level)
            {
                break;
            }
            if constexpr (RecordStatsLevel > 1)
            {
                m_metrics.candidates_inspected++;
            }
            if (!m_active[i])
            {
                continue;
            }
            const auto& [start, end] = m_records[i].interval;
            // check if they intersect
            if (intervals_touch_or_cross(query, {start, end}))
            {
                if constexpr (RecordStatsLevel > 0)
                {
                    m_metrics.callbacks++;
                }
                callback(i, {start, end});
            }
        }
    }

    void popInterval(Interval interval) override
    {
        for (Index index = 0; index < m_records.size(); ++index)
        {
            if (m_active[index] && m_records[index].interval == interval)
            {
                popInterval(index);
                return;
            }
        }
    }
    void popInterval(Index idx) override
    {
        if (idx < m_active.size() && m_active[idx])
        {
            if constexpr (RecordStatsLevel > 0)
            {
                m_metrics.intervals_popped++;
            }
            m_active[idx] = false;
        }
    }

    void activateInterval(Index idx, size_t level) override
    {
        if (idx < m_active.size())
        {
            m_records[idx].level = level;
            m_active[idx] = true;
        }
    }

    void clear() override
    {
        std::fill(m_active.begin(), m_active.end(), false);
    }

    void reset() override
    {
        std::fill(m_active.begin(), m_active.end(), true);
        for (Index index = 0; index < m_records.size(); ++index)
        {
            m_records[index].level = m_initial_levels[index];
        }
    }

    IntersectionIndexMetrics emit_metrics() const override
    {
        return m_metrics;
    }

private:
    std::vector<IntersectionRecord> m_records;
    std::vector<char> m_active;
    std::vector<size_t> m_initial_levels;
    bool m_incremental{false};
    mutable IntersectionIndexMetrics m_metrics;
};
