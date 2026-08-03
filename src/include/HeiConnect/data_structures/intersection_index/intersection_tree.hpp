#pragma once

#include "HeiConnect/data_structures/intersection_index/baseline.hpp"

#include <algorithm>
#include <limits>
#include <numeric>

template<int RecordStatsLevel = 0>
class IntersectionTreeIdx : public BaseIntersectionIdx<RecordStatsLevel>
{
public:
    using typename BaseIntersectionIdx<RecordStatsLevel>::Index;
    using typename BaseIntersectionIdx<RecordStatsLevel>::Interval;

    IntersectionTreeIdx() = default;

    explicit IntersectionTreeIdx(const std::vector<Interval>& intervals) : m_intervals(intervals)
    {
        m_active.resize(intervals.size(), true);
        m_node_by_interval.resize(intervals.size());
        std::vector<Index> indices(intervals.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](Index left, Index right) {
            return m_intervals[left] < m_intervals[right];
        });
        m_root = build(indices, 0, indices.size(), -1);
        m_initial_max_end.reserve(m_nodes.size());
        for (const Node& node : m_nodes)
        {
            m_initial_max_end.push_back(node.max_end);
        }
    }

    std::unique_ptr<BaseIntersectionIdx<RecordStatsLevel>> make(const std::vector<Interval>& intervals) const override
    {
        return std::make_unique<IntersectionTreeIdx<RecordStatsLevel>>(intervals);
    }

    void forEachIntersection(std::function<void(Index, Interval)> callback, Interval query) const override
    {
        if constexpr (RecordStatsLevel > 0)
        {
            m_metrics.queries++;
        }
        queryTree(m_root, query, callback);
    }

    void popInterval(Interval interval) override
    {
        for (Index index = 0; index < m_intervals.size(); ++index)
        {
            if (m_active[index] && m_intervals[index] == interval)
            {
                popInterval(index);
                return;
            }
        }
    }

    void popInterval(Index index) override
    {
        if (index < m_active.size() && m_active[index])
        {
            if constexpr (RecordStatsLevel > 0)
            {
                m_metrics.intervals_popped++;
            }
            m_active[index] = false;
            int node = m_node_by_interval[index];
            while (node >= 0)
            {
                update(node);
                node = m_nodes[node].parent;
            }
        }
    }

    void reset() override
    {
        std::fill(m_active.begin(), m_active.end(), true);
        for (size_t node = 0; node < m_nodes.size(); ++node)
        {
            m_nodes[node].max_end = m_initial_max_end[node];
        }
    }

    IntersectionIndexMetrics emit_metrics() const override
    {
        return m_metrics;
    }

private:
    struct Node
    {
        Index interval;
        size_t max_end;
        int left{-1};
        int right{-1};
        int parent{-1};
    };

    static constexpr size_t inactive_end = std::numeric_limits<size_t>::min();

    int build(const std::vector<Index>& indices, size_t begin, size_t end, int parent)
    {
        if (begin == end)
        {
            return -1;
        }

        const size_t middle = begin + (end - begin) / 2;
        const int node = static_cast<int>(m_nodes.size());
        m_nodes.push_back({indices[middle], std::get<1>(m_intervals[indices[middle]]), -1, -1, parent});
        m_node_by_interval[indices[middle]] = node;
        m_nodes[node].left = build(indices, begin, middle, node);
        m_nodes[node].right = build(indices, middle + 1, end, node);
        update(node);
        return node;
    }

    size_t maxEnd(int node) const
    {
        return node < 0 ? inactive_end : m_nodes[node].max_end;
    }

    void update(int node)
    {
        const auto index = m_nodes[node].interval;
        const size_t own_end = m_active[index] ? std::get<1>(m_intervals[index]) : inactive_end;
        m_nodes[node].max_end = std::max({own_end, maxEnd(m_nodes[node].left), maxEnd(m_nodes[node].right)});
    }

    static bool intersects(Interval left, Interval right)
    {
        const auto [a, b] = left;
        const auto [c, d] = right;
        return a == c || a == d || b == c || b == d || (a < c && c < b && b < d) || (c < a && a < d && d < b);
    }

    void queryTree(int node, Interval query, const std::function<void(Index, Interval)>& callback) const
    {
        if (node < 0 || m_nodes[node].max_end < std::get<0>(query))
        {
            return;
        }

        if constexpr (RecordStatsLevel > 1)
        {
            m_metrics.candidates_inspected++;
        }

        queryTree(m_nodes[node].left, query, callback);
        const Index index = m_nodes[node].interval;
        const Interval interval = m_intervals[index];
        if (m_active[index] && std::get<0>(interval) <= std::get<1>(query) && intersects(interval, query))
        {
            if constexpr (RecordStatsLevel > 0)
            {
                m_metrics.callbacks++;
            }
            callback(index, interval);
        }
        if (std::get<0>(interval) <= std::get<1>(query))
        {
            queryTree(m_nodes[node].right, query, callback);
        }
    }

    std::vector<Interval> m_intervals;
    std::vector<Node> m_nodes;
    std::vector<char> m_active;
    std::vector<int> m_node_by_interval;
    std::vector<size_t> m_initial_max_end;
    int m_root{-1};
    mutable IntersectionIndexMetrics m_metrics;
};
