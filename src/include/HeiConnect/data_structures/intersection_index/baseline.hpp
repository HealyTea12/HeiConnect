#pragma once

#include <algorithm>
#include <functional>
#include <memory>
#include <vector>
#include <tuple>

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
    std::tuple<size_t, size_t> interval;
    size_t level;
};

template<int RecordStatsLevel = 0>
class BaseIntersectionIdx
{
public:
    virtual ~BaseIntersectionIdx() = default;
    using Interval = std::tuple<size_t, size_t>;
    using Index = size_t;

    virtual std::unique_ptr<BaseIntersectionIdx<RecordStatsLevel>> make(
        const std::vector<IntersectionRecord>& records) const = 0;

    virtual void forEachIntersection(
        std::function<void(Index, Interval)> callback,
        Interval query,
        size_t exclusive_level) const = 0;

    virtual void popInterval(Interval interval) = 0;
    virtual void popInterval(Index idx) = 0;
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
    }

    std::unique_ptr<BaseIntersectionIdx<RecordStatsLevel>> make(
        const std::vector<IntersectionRecord>& records) const override
    {
        return std::make_unique<BaselineIntersectionIdx<RecordStatsLevel>>(records);
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
        auto [q_start, q_end] = query;
        for (size_t i{0}; i < m_records.size(); ++i)
        {
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
            auto intersection_size = std::min(q_end, end) - std::max(q_start, start);
            auto query_size = q_end - q_start;
            auto interval_size = end - start;
            bool touch_endpoints = (q_start == end) || (q_start == start) || (q_end == end) || (q_end == start);
            if (intersection_size < std::min(query_size, interval_size) || touch_endpoints)
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

    void reset() override
    {
        std::fill(m_active.begin(), m_active.end(), true);
    }

    IntersectionIndexMetrics emit_metrics() const override
    {
        return m_metrics;
    }

private:
    std::vector<IntersectionRecord> m_records;
    std::vector<char> m_active;
    mutable IntersectionIndexMetrics m_metrics;
};
