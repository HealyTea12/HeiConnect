#pragma once

#include <functional>
#include <memory>
#include <vector>
#include <tuple>

class BaseIntersectionIdx
{
public:
    virtual ~BaseIntersectionIdx() = default;
    using Interval = std::tuple<size_t, size_t>;
    using Index = size_t;

    virtual void forEachIntersection(std::function<void(Index, Interval)> callback, Interval query) const = 0;

    virtual void popInterval(Interval interval) = 0;
    virtual void popInterval(Index idx) = 0;
};

class BaselineIntersectionIdx : public BaseIntersectionIdx
{
public:
    BaselineIntersectionIdx(const std::vector<std::tuple<size_t, size_t>>& intervals)
    {
        m_intervals = std::make_shared<std::vector<std::tuple<size_t, size_t>>>(intervals);
    }

    void forEachIntersection(std::function<void(Index, Interval)> callback, Interval query) const override
    {
        auto [q_start, q_end] = query;
        for (size_t i{0}; i < m_intervals->size(); ++i)
        {
            const auto& [start, end] = (*m_intervals)[i];
            // check if they intersect
            auto intersection_size = std::min(q_end, end) - std::max(q_start, start);
            auto query_size = q_end - q_start;
            auto interval_size = end - start;
            bool touch_endpoints = (q_start == end) || (q_start == start) || (q_end == end) || (q_end == start);
            if (intersection_size < std::min(query_size, interval_size) || touch_endpoints)
            {
                callback(i, {start, end});
            }
        }
    }

    void popInterval(Interval interval) override
    {
        return;
    }
    void popInterval(Index idx) override
    {
        return;
    }

private:
    std::shared_ptr<std::vector<std::tuple<size_t, size_t>>> m_intervals;
};