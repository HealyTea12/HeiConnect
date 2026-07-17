#pragma once

#include "HeiConnect/data_structures/intersection_index/baseline.hpp"

#include <stdexcept>
#include <algorithm>
#include <tuple>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <utility>
#include <vector>


/**
 * Cycle Domination Baseline.
 * Inputs:
 *  - links: a vector of tuples of [u,v,w] representing unique links on a cycle graph.
 *           They are normalized so that 0 <= u < v < n_nodes.
 * - n_nodes: the number of nodes in the cycle graph.
 *
 * This algorithm will performa a "bounded" Dijkstra from each node in the cycle graph.
 * The goal of the Dijkstra is to calculate the shortest distance to each other node in the cycle, following the link
 * graph. The link intersection graph is a graph where conceptually each node represents a link where you pay its weight
 * to enter it. There are also nodes representing the cycle nodes, which have 0 cost.
 *
 * Hence, this algorithm performs a Dijkstra on the link intersection graph starting at each cycle node.
 * It calculates the maximum weight of a link incident to the source node and uses it as a cutoff.
 * Once it relaxes a link, it may pop all neighbouring links off the intersectionIndex, and their final distance is
 * calculated. Once it relaxes a link, it has finalized the distance to that link's endpoints, so it can mark them as
 * complete and their distance is calculated. Once all endpoints are complete, the Dijkstra can terminate early.
 *
 * After each run of Dijkstra, we check if any of the links incident to the source node have a shorter distance to their
 * other endpoint. if so, we mark that link as removable.
 */
template<typename NodeID, typename Weight, typename IntersectionIdx = BaselineIntersectionIdx>
auto cycle_domination_baseline(const std::vector<std::tuple<NodeID, NodeID, Weight>>& links, uint64_t n_nodes)
{
    using LinkID = size_t;
    using QueueEntry = std::pair<Weight, LinkID>;

    const LinkID invalid_link = std::numeric_limits<LinkID>::max();

    std::vector<std::vector<LinkID>> incident(n_nodes);
    std::vector<std::tuple<size_t, size_t>> intervals;
    intervals.reserve(links.size());

    for (LinkID id = 0; id < links.size(); ++id)
    {
        const auto& [u, v, weight] = links[id];
        incident[static_cast<size_t>(u)].push_back(id);
        incident[static_cast<size_t>(v)].push_back(id);
        intervals.emplace_back(static_cast<size_t>(u), static_cast<size_t>(v));
    }

    std::vector<char> removable(links.size(), false);
    std::vector<char> explored(links.size());
    std::vector<char> complete(n_nodes);
    std::vector<Weight> vertex_distance(n_nodes);
    std::vector<Weight> value_to_beat(n_nodes);
    std::vector<LinkID> candidate_link(n_nodes);

    for (size_t source = 0; source < n_nodes; ++source)
    {
        if (incident[source].empty())
        {
            continue;
        }

        std::fill(explored.begin(), explored.end(), false);
        std::fill(complete.begin(), complete.end(), false);
        std::fill(candidate_link.begin(), candidate_link.end(), invalid_link);
        size_t n_complete_vertices = 0;

        Weight cutoff{};
        for (LinkID id : incident[source])
        {
            const auto& [u, v, weight] = links[id];
            const size_t target = static_cast<size_t>(u) == source ? static_cast<size_t>(v) : static_cast<size_t>(u);
            value_to_beat[target] = weight;
            candidate_link[target] = id;
            cutoff = std::max(cutoff, weight);
        }

        IntersectionIdx intersection_index(intervals);
        std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;

        auto complete_endpoint = [&](NodeID endpoint, Weight value) {
            const size_t vertex = static_cast<size_t>(endpoint);
            if (complete[vertex])
            {
                return;
            }

            n_complete_vertices++;
            complete[vertex] = true;
            vertex_distance[vertex] = value;

            if (candidate_link[vertex] != invalid_link && vertex_distance[vertex] < value_to_beat[vertex])
            {
                removable[candidate_link[vertex]] = true;
            }
        };

        for (LinkID id : incident[source])
        {
            explored[id] = true;
            intersection_index.popInterval(id);
            queue.emplace(std::get<2>(links[id]), id);
        }

        while (!queue.empty())
        {
            const auto [current_distance, current] = queue.top();
            queue.pop();

            const auto& [u, v, weight] = links[current];
            complete_endpoint(u, current_distance);
            complete_endpoint(v, current_distance);

            if (current_distance >= cutoff || n_complete_vertices == n_nodes)
            {
                break;
            }

            std::vector<LinkID> neighbours;
            intersection_index.forEachIntersection(
                [&](LinkID next, typename IntersectionIdx::Interval) {
                    const Weight next_weight = std::get<2>(links[next]);
                    if (!explored[next] && next_weight < cutoff - current_distance)
                    {
                        explored[next] = true;
                        neighbours.push_back(next);
                    }
                },
                intervals[current]);

            for (LinkID next : neighbours)
            {
                intersection_index.popInterval(next);
                queue.emplace(current_distance + std::get<2>(links[next]), next);
            }
        }
    }

    std::vector<int> result;
    for (LinkID id = 0; id < links.size(); ++id)
    {
        if (removable[id])
        {
            result.push_back(static_cast<int>(id));
        }
    }

    return result;
}

// /*
// template<typename NodeID, typename Weight>
// auto full_cycle_link_domination(const std::vector<std::tuple<NodeID, NodeID, Weight>>& links, uint64_t n_nodes)
// {
//     auto comp = [](const auto& a, const auto& b) {
//         return std::get<0>(a) < std::get<0>(b);
//     };
//     std::priority_queue<std::tuple<Weight, NodeID, NodeID>> link_queue;
//     std::sort(links.begin(), links.end(), comp);
//     for (uint64_t i{}; i < n_nodes; i++)
//     {
//         Weight upper_limit = 0;
//         auto start = std::lower_bound(links.begin(), links.end(), std::make_tuple(i, NodeID{}, Weight{}), comp);
//         auto end = std::upper_bound(links.begin(), links.end(), std::make_tuple(i, NodeID{}, Weight{}), comp);
//         for (auto it = start; it != end; it++)
//         {
//             const auto& [u, v, w] = *it;
//             if (w > upper_limit)
//             {
//                 upper_limit = w;
//             }
//         }
//     }
//     return;
// }
//
// template<typename NodeID, typename Weight>
// struct Link
// {
//     NodeID u;
//     NodeID v;
//     Weight weight;
// };
//
// /**
//  * Returns the input indices of all links l = {a,b} for which there is
//  * an a-to-b chain of touching links with total weight *strictly* less
//  * than weight(l).
//  *
//  * Two links touch if they share an endpoint or cross on the cycle.
//  * i.e l1 = {a,b} and l2 = {c,d} touch if a=c, a=d, b=c, b=d, or a<c<b<d or c<a<d<b.
//  *
//  * Assumptions:
//  *   - 0 <= link.u < link.v < n;
//  *   - links are unique;
//  *   - weights are nonnegative.
//  */
// template<typename NodeID, typename Weight>
// std::vector<int> cycle_link_domination_simple_baseline(int n, const std::vector<Link<NodeID, Weight>>& links)
// {
//     constexpr Weight INF = std::numeric_limits<Weight>::max();
//
//     constexpr Weight NO_TARGET = -1;
//     const size_t INVALID_LINK_ID = std::numeric_limits<size_t>::max();
//
//     const size_t m = links.size();
//
//     /*
//      * Link ids of links incident to each cycle vertex.
//      * TODO: Flatten into csr
//      */
//     std::vector<std::vector<size_t>> incident(n);
//
//     for (size_t id = 0; id < m; ++id)
//     {
//         incident[links[id].u].push_back(id);
//         incident[links[id].v].push_back(id);
//     }
//
//     /*
//      * Link indices sorted by their smaller endpoint.
//      */
//     std::vector<size_t> sortedLinks(m);
//     std::iota(sortedLinks.begin(), sortedLinks.end(), 0ull);
//
//     std::sort(sortedLinks.begin(), sortedLinks.end(), [&](int first, int second) {
//         if (links[first].u != links[second].u)
//         {
//             return links[first].u < links[second].u;
//         }
//
//         return links[first].v < links[second].v;
//     });
//
//     auto touch = [&](size_t first, size_t second) {
//         const CycleLink& a = links[first];
//         const CycleLink& b = links[second];
//
//         const bool shareEndpoint = a.u == b.u || a.u == b.v || a.v == b.u || a.v == b.v;
//
//         const bool cross = (a.u < b.u && b.u < a.v && a.v < b.v) || (b.u < a.u && a.u < b.v && b.v < a.v);
//
//         return shareEndpoint || cross;
//     };
//
//     std::vector<char> removable(m, false);
//
//     std::vector<Weight> linkDistance(m);
//
//     /*
//      * For a fixed source s:
//      *
//      * valueToBeat[v] = weight of the unique link {s,v},
//      *                  or NO_TARGET if there is no such link.
//      *
//      * candidateLink[v] = input index of {s,v}.
//      */
//     std::vector<Weight> valueToBeat(n);
//     std::vector<int> candidateLink(n);
//
//     using QueueEntry = std::pair<Weight, size_t>; // (distance, link id)
//
//     /*
//      * Run the link-graph Dijkstra once from every cycle vertex that is
//      * incident to at least one link.
//      */
//     for (size_t source = 0; source < n; ++source)
//     {
//         if (incident[source].empty())
//         {
//             continue;
//         }
//         // std::fill(linkDistance.begin(), linkDistance.end(), INF);
//         std::fill(valueToBeat.begin(), valueToBeat.end(), NO_TARGET);
//         std::fill(candidateLink.begin(), candidateLink.end(), INVALID_LINK_ID);
//         /*
//          * No distance greater than or equal to cutoff can eliminate
//          * any link incident to this source.
//          */
//         Weight cutoff = 0;
//
//         for (size_t id : incident[source])
//         {
//             // which one is target/source setting both might be faster, but who cares
//             const size_t target;
//             if (links[id].u == source)
//             {
//                 target = links[id].v;
//             }
//             else
//             {
//                 target = links[id].u;
//             }
//
//             valueToBeat[target] = links[id].weight;
//             candidateLink[target] = id;
//
//             cutoff = std::max(cutoff, links[id].weight);
//         }
//
//         std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> pq;
//
//         /*
//          * Reaching a link also reaches both of its endpoints.
//          *
//          * The distance need not be finalized before this check: every
//          * tentative label corresponds to an actual link chain, so any
//          * value strictly below valueToBeat[vertex] is already a valid
//          * certificate of removability.
//          * After relaxing an endpoint once, it is actually never necessary to relax it again.
//          * Once all endpoints have been relaxed, the algorithm can terminate early.
//          */
//         auto relaxEndpoint = [&](NodeID vertex, Weight value) {
//             if (valueToBeat[vertex] != NO_TARGET && value < valueToBeat[vertex])
//             {
//                 removable[candidateLink[vertex]] = true;
//             }
//         };
//
//         auto relaxLink = [&](size_t id, Weight value) {
//             /*
//              * Values at least cutoff cannot strictly beat any candidate
//              * link belonging to this source run.
//              */
//             if (value >= cutoff || value >= linkDistance[id])
//             {
//                 return;
//             }
//
//             linkDistance[id] = value;
//             pq.emplace(value, id);
//
//             relaxEndpoint(links[id].u, value);
//             relaxEndpoint(links[id].v, value);
//         };
//
//         /*
//          * The virtual source vertex reaches an incident link by paying
//          * that link's weight.
//          */
//         for (int id : incident[source])
//         {
//             relaxLink(id, links[id].weight);
//         }
//
//         /*
//          * Continue only while the minimum queue value can still beat at
//          * least one direct link incident to the source.
//          */
//         while (!pq.empty() && pq.top().first < cutoff)
//         {
//             const auto [distance, current] = pq.top();
//             pq.pop();
//
//             if (distance != linkDistance[current])
//             {
//                 continue; // Stale queue entry.
//             }
//
//             const CycleLink& currentLink = links[current];
//
//             /*
//              * Examine potential neighbors in increasing order of their
//              * smaller endpoint.
//              *
//              * When next.u > currentLink.v, next cannot cross currentLink
//              * and cannot share an endpoint with it, so the scan can stop.
//              */
//             for (int next : sortedLinks)
//             {
//                 const CycleLink& nextLink = links[next];
//
//                 if (nextLink.u > currentLink.v)
//                 {
//                     break;
//                 }
//
//                 if (next == current || !touch(current, next))
//                 {
//                     continue;
//                 }
//
//                 /*
//                  * Require
//                  *
//                  *     distance + nextLink.weight < cutoff.
//                  *
//                  * This form avoids integer overflow.
//                  */
//                 if (nextLink.weight >= cutoff - distance)
//                 {
//                     continue;
//                 }
//
//                 relaxLink(next, distance + nextLink.weight);
//             }
//         }
//     }
//
//     std::vector<int> result;
//
//     for (int id = 0; id < m; ++id)
//     {
//         if (removable[id])
//         {
//             result.push_back(id);
//         }
//     }
//
//     return result;
// }
//
//
// struct CycleLink
// {
//     int u;
//     int v;
//     std::int64_t weight;
// };
//
// /**
//  * Cycle link reduction algorithm based on bounded dijkstra + segment trees for calculating neighbours on
//  * the cycle intersection graph.
//  */
// namespace cycle_link_reduction_detail
// {
//
//     class MaxReporter
//     {
//     public:
//         explicit MaxReporter(int size) : size_(size), base_(1)
//         {
//             while (base_ < std::max(1, size_))
//             {
//                 base_ <<= 1;
//             }
//             tree_.assign(2 * base_, kNegativeInfinity);
//         }
//
//         // Leaves are ordered by the left endpoint of each link.
//         // A leaf stores the right endpoint of the link if it is active.
//         void build(const std::vector<int>& order, const std::vector<char>& active, const std::vector<int>& right)
//         {
//             std::fill(tree_.begin(), tree_.end(), kNegativeInfinity);
//
//             for (int pos = 0; pos < size_; ++pos)
//             {
//                 const int id = order[pos];
//                 if (active[id])
//                 {
//                     tree_[base_ + pos] = right[id];
//                 }
//             }
//
//             for (int node = base_ - 1; node >= 1; --node)
//             {
//                 tree_[node] = std::max(tree_[2 * node], tree_[2 * node + 1]);
//             }
//         }
//
//         void eraseAt(int position)
//         {
//             int node = base_ + position;
//
//             if (tree_[node] == kNegativeInfinity)
//             {
//                 return;
//             }
//
//             tree_[node] = kNegativeInfinity;
//
//             for (node >>= 1; node >= 1; node >>= 1)
//             {
//                 const int newValue = std::max(tree_[2 * node], tree_[2 * node + 1]);
//
//                 if (tree_[node] == newValue)
//                 {
//                     break;
//                 }
//
//                 tree_[node] = newValue;
//             }
//         }
//
//         // Returns a position in [queryLeft, queryRight) whose stored value is
//         // strictly larger than threshold, or -1 if none exists.
//         int findPosition(int queryLeft, int queryRight, int threshold) const
//         {
//             if (queryLeft >= queryRight || tree_[1] <= threshold)
//             {
//                 return -1;
//             }
//
//             return findPositionRecursive(1, 0, base_, queryLeft, queryRight, threshold);
//         }
//
//     private:
//         static constexpr int kNegativeInfinity = -1;
//
//         int size_;
//         int base_;
//         std::vector<int> tree_;
//
//         int
//         findPositionRecursive(int node, int nodeLeft, int nodeRight, int queryLeft, int queryRight, int threshold)
//         const
//         {
//             if (nodeRight <= queryLeft || queryRight <= nodeLeft || tree_[node] <= threshold)
//             {
//                 return -1;
//             }
//
//             if (nodeRight - nodeLeft == 1)
//             {
//                 return nodeLeft < size_ ? nodeLeft : -1;
//             }
//
//             const int middle = nodeLeft + (nodeRight - nodeLeft) / 2;
//
//             int answer = findPositionRecursive(2 * node, nodeLeft, middle, queryLeft, queryRight, threshold);
//
//             if (answer != -1)
//             {
//                 return answer;
//             }
//
//             return findPositionRecursive(2 * node + 1, middle, nodeRight, queryLeft, queryRight, threshold);
//         }
//     };
//
//     class MinReporter
//     {
//     public:
//         explicit MinReporter(int size) : size_(size), base_(1)
//         {
//             while (base_ < std::max(1, size_))
//             {
//                 base_ <<= 1;
//             }
//             tree_.assign(2 * base_, kInfinity);
//         }
//
//         // Leaves are ordered by the right endpoint of each link.
//         // A leaf stores the left endpoint of the link if it is active.
//         void build(const std::vector<int>& order, const std::vector<char>& active, const std::vector<int>& left)
//         {
//             std::fill(tree_.begin(), tree_.end(), kInfinity);
//
//             for (int pos = 0; pos < size_; ++pos)
//             {
//                 const int id = order[pos];
//                 if (active[id])
//                 {
//                     tree_[base_ + pos] = left[id];
//                 }
//             }
//
//             for (int node = base_ - 1; node >= 1; --node)
//             {
//                 tree_[node] = std::min(tree_[2 * node], tree_[2 * node + 1]);
//             }
//         }
//
//         void eraseAt(int position)
//         {
//             int node = base_ + position;
//
//             if (tree_[node] == kInfinity)
//             {
//                 return;
//             }
//
//             tree_[node] = kInfinity;
//
//             for (node >>= 1; node >= 1; node >>= 1)
//             {
//                 const int newValue = std::min(tree_[2 * node], tree_[2 * node + 1]);
//
//                 if (tree_[node] == newValue)
//                 {
//                     break;
//                 }
//
//                 tree_[node] = newValue;
//             }
//         }
//
//         // Returns a position in [queryLeft, queryRight) whose stored value is
//         // strictly smaller than threshold, or -1 if none exists.
//         int findPosition(int queryLeft, int queryRight, int threshold) const
//         {
//             if (queryLeft >= queryRight || tree_[1] >= threshold)
//             {
//                 return -1;
//             }
//
//             return findPositionRecursive(1, 0, base_, queryLeft, queryRight, threshold);
//         }
//
//     private:
//         static constexpr int kInfinity = std::numeric_limits<int>::max();
//
//         int size_;
//         int base_;
//         std::vector<int> tree_;
//
//         int
//         findPositionRecursive(int node, int nodeLeft, int nodeRight, int queryLeft, int queryRight, int threshold)
//         const
//         {
//             if (nodeRight <= queryLeft || queryRight <= nodeLeft || tree_[node] >= threshold)
//             {
//                 return -1;
//             }
//
//             if (nodeRight - nodeLeft == 1)
//             {
//                 return nodeLeft < size_ ? nodeLeft : -1;
//             }
//
//             const int middle = nodeLeft + (nodeRight - nodeLeft) / 2;
//
//             int answer = findPositionRecursive(2 * node, nodeLeft, middle, queryLeft, queryRight, threshold);
//
//             if (answer != -1)
//             {
//                 return answer;
//             }
//
//             return findPositionRecursive(2 * node + 1, middle, nodeRight, queryLeft, queryRight, threshold);
//         }
//     };
//
// } // namespace cycle_link_reduction_detail
//
// /**
//  * Returns the indices of all strictly removable links.
//  *
//  * A link l = {a,b} is returned exactly when there is an a-to-b chain
//  * of links having total weight strictly smaller than l.weight.
//  *
//  * Consecutive links in a chain must either:
//  *   - cross as chords of the cycle, or
//  *   - share an endpoint.
//  *
//  * Requirements:
//  *   - vertices are numbered 0,...,n-1 in cyclic order;
//  *   - link weights are nonnegative;
//  *   - links have distinct endpoints.
//  *
//  * The returned indices are in increasing input order.
//  */
// std::vector<int> removableCycleLinks(int n, const std::vector<CycleLink>& links)
// {
//     using namespace cycle_link_reduction_detail;
//     using Distance = std::int64_t;
//
//     constexpr Distance kInfinity = std::numeric_limits<Distance>::max();
//
//     const int numberOfLinks = static_cast<int>(links.size());
//
//     if (n < 3)
//     {
//         throw std::invalid_argument("A cycle must contain at least three vertices.");
//     }
//
//     if (numberOfLinks == 0)
//     {
//         return {};
//     }
//
//     /*
//      * Coordinate-compress only vertices appearing as link endpoints.
//      * Thus the storage is O(numberOfLinks), even when n is very large.
//      */
//     std::vector<int> coordinates;
//     coordinates.reserve(2 * numberOfLinks);
//
//     for (const CycleLink& link : links)
//     {
//         if (link.u < 0 || link.u >= n || link.v < 0 || link.v >= n)
//         {
//             throw std::invalid_argument("A link endpoint is outside [0,n).");
//         }
//
//         if (link.u == link.v)
//         {
//             throw std::invalid_argument("Self-links are not supported.");
//         }
//
//         if (link.weight < 0)
//         {
//             throw std::invalid_argument("Link weights must be nonnegative.");
//         }
//
//         coordinates.push_back(link.u);
//         coordinates.push_back(link.v);
//     }
//
//     std::sort(coordinates.begin(), coordinates.end());
//     coordinates.erase(std::unique(coordinates.begin(), coordinates.end()), coordinates.end());
//
//     const int numberOfEndpointVertices = static_cast<int>(coordinates.size());
//
//     std::vector<int> left(numberOfLinks);
//     std::vector<int> right(numberOfLinks);
//     std::vector<Distance> weight(numberOfLinks);
//     std::vector<int> degree(numberOfEndpointVertices, 0);
//
//     for (int id = 0; id < numberOfLinks; ++id)
//     {
//         int a = static_cast<int>(
//             std::lower_bound(coordinates.begin(), coordinates.end(), links[id].u) - coordinates.begin());
//
//         int b = static_cast<int>(
//             std::lower_bound(coordinates.begin(), coordinates.end(), links[id].v) - coordinates.begin());
//
//         if (a > b)
//         {
//             std::swap(a, b);
//         }
//
//         left[id] = a;
//         right[id] = b;
//         weight[id] = links[id].weight;
//
//         ++degree[a];
//         ++degree[b];
//     }
//
//     /*
//      * Store endpoint-to-link incidence lists in CSR form.
//      * This represents shared-endpoint adjacencies without materializing
//      * the corresponding cliques.
//      */
//     std::vector<int> incidenceStart(numberOfEndpointVertices + 1, 0);
//
//     for (int vertex = 0; vertex < numberOfEndpointVertices; ++vertex)
//     {
//         incidenceStart[vertex + 1] = incidenceStart[vertex] + degree[vertex];
//     }
//
//     std::vector<int> incidence(2 * numberOfLinks);
//     std::vector<int> incidenceCursor = incidenceStart;
//
//     for (int id = 0; id < numberOfLinks; ++id)
//     {
//         incidence[incidenceCursor[left[id]]++] = id;
//         incidence[incidenceCursor[right[id]]++] = id;
//     }
//
//     /*
//      * The two orderings support implicit crossing-neighbor reporting.
//      */
//     std::vector<int> byLeft(numberOfLinks);
//     std::vector<int> byRight(numberOfLinks);
//     std::vector<int> positionByLeft(numberOfLinks);
//     std::vector<int> positionByRight(numberOfLinks);
//
//     for (int id = 0; id < numberOfLinks; ++id)
//     {
//         byLeft[id] = id;
//         byRight[id] = id;
//     }
//
//     std::sort(byLeft.begin(), byLeft.end(), [&](int first, int second) {
//         return std::tie(left[first], right[first], first) < std::tie(left[second], right[second], second);
//     });
//
//     std::sort(byRight.begin(), byRight.end(), [&](int first, int second) {
//         return std::tie(right[first], left[first], first) < std::tie(right[second], left[second], second);
//     });
//
//     for (int position = 0; position < numberOfLinks; ++position)
//     {
//         positionByLeft[byLeft[position]] = position;
//         positionByRight[byRight[position]] = position;
//     }
//
//     // First position in byLeft whose left endpoint is > value.
//     auto upperLeft = [&](int value) {
//         int low = 0;
//         int high = numberOfLinks;
//
//         while (low < high)
//         {
//             const int middle = low + (high - low) / 2;
//
//             if (left[byLeft[middle]] <= value)
//             {
//                 low = middle + 1;
//             }
//             else
//             {
//                 high = middle;
//             }
//         }
//
//         return low;
//     };
//
//     // First position in byLeft whose left endpoint is >= value.
//     auto lowerLeft = [&](int value) {
//         int low = 0;
//         int high = numberOfLinks;
//
//         while (low < high)
//         {
//             const int middle = low + (high - low) / 2;
//
//             if (left[byLeft[middle]] < value)
//             {
//                 low = middle + 1;
//             }
//             else
//             {
//                 high = middle;
//             }
//         }
//
//         return low;
//     };
//
//     // First position in byRight whose right endpoint is > value.
//     auto upperRight = [&](int value) {
//         int low = 0;
//         int high = numberOfLinks;
//
//         while (low < high)
//         {
//             const int middle = low + (high - low) / 2;
//
//             if (right[byRight[middle]] <= value)
//             {
//                 low = middle + 1;
//             }
//             else
//             {
//                 high = middle;
//             }
//         }
//
//         return low;
//     };
//
//     // First position in byRight whose right endpoint is >= value.
//     auto lowerRight = [&](int value) {
//         int low = 0;
//         int high = numberOfLinks;
//
//         while (low < high)
//         {
//             const int middle = low + (high - low) / 2;
//
//             if (right[byRight[middle]] < value)
//             {
//                 low = middle + 1;
//             }
//             else
//             {
//                 high = middle;
//             }
//         }
//
//         return low;
//     };
//
//     /*
//      * One shortest-path run from a vertex tests every link assigned to
//      * that vertex. Assigning a link to its higher-degree endpoint is a
//      * heuristic that reduces the number of source runs on instances
//      * such as stars.
//      */
//     auto sourceOf = [&](int id) {
//         if (degree[left[id]] != degree[right[id]])
//         {
//             return degree[left[id]] > degree[right[id]] ? left[id] : right[id];
//         }
//
//         return left[id];
//     };
//
//     auto targetOf = [&](int id, int source) {
//         return left[id] == source ? right[id] : left[id];
//     };
//
//     std::vector<int> testOrder(numberOfLinks);
//
//     for (int id = 0; id < numberOfLinks; ++id)
//     {
//         testOrder[id] = id;
//     }
//
//     std::sort(testOrder.begin(), testOrder.end(), [&](int first, int second) {
//         const int firstSource = sourceOf(first);
//         const int secondSource = sourceOf(second);
//
//         if (firstSource != secondSource)
//         {
//             return firstSource < secondSource;
//         }
//
//         return first < second;
//     });
//
//     MaxReporter maxRightReporter(numberOfLinks);
//     MinReporter minLeftReporter(numberOfLinks);
//
//     std::vector<char> active(numberOfLinks, 0);
//     std::vector<Distance> endpointDistance(numberOfEndpointVertices, kInfinity);
//     std::vector<char> removable(numberOfLinks, 0);
//
//     int groupBegin = 0;
//
//     while (groupBegin < numberOfLinks)
//     {
//         const int source = sourceOf(testOrder[groupBegin]);
//
//         int groupEnd = groupBegin;
//         Distance cutoff = 0;
//
//         while (groupEnd < numberOfLinks && sourceOf(testOrder[groupEnd]) == source)
//         {
//             cutoff = std::max(cutoff, weight[testOrder[groupEnd]]);
//             ++groupEnd;
//         }
//
//         /*
//          * We only need paths with total weight < cutoff.
//          * Therefore links of weight >= cutoff cannot be part of a
//          * relevant path.
//          */
//         if (cutoff > 0)
//         {
//             for (int id = 0; id < numberOfLinks; ++id)
//             {
//                 active[id] = static_cast<char>(weight[id] < cutoff);
//             }
//
//             maxRightReporter.build(byLeft, active, right);
//
//             minLeftReporter.build(byRight, active, left);
//
//             std::fill(endpointDistance.begin(), endpointDistance.end(), kInfinity);
//
//             using QueueItem = std::pair<Distance, int>;
//
//             std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<QueueItem>> queue;
//
//             /*
//              * The first settled neighbor that discovers a link gives
//              * its final shortest-path label: every transition entering
//              * link id adds exactly weight[id].
//              */
//             auto discover = [&](int id, Distance predecessorDistance) {
//                 if (!active[id])
//                 {
//                     return;
//                 }
//
//                 active[id] = 0;
//
//                 maxRightReporter.eraseAt(positionByLeft[id]);
//
//                 minLeftReporter.eraseAt(positionByRight[id]);
//
//                 /*
//                  * Only distances strictly below cutoff matter.
//                  * Writing the test this way also avoids overflow.
//                  */
//                 if (weight[id] < cutoff - predecessorDistance)
//                 {
//                     queue.emplace(predecessorDistance + weight[id], id);
//                 }
//             };
//
//             /*
//              * Expanding an endpoint implicitly traverses all
//              * shared-endpoint adjacencies.
//              */
//             auto expandEndpoint = [&](int vertex, Distance distance) {
//                 if (endpointDistance[vertex] != kInfinity)
//                 {
//                     return;
//                 }
//
//                 endpointDistance[vertex] = distance;
//
//                 for (int position = incidenceStart[vertex]; position < incidenceStart[vertex + 1]; ++position)
//                 {
//                     discover(incidence[position], distance);
//                 }
//             };
//
//             expandEndpoint(source, 0);
//
//             while (!queue.empty())
//             {
//                 const auto [distance, current] = queue.top();
//                 queue.pop();
//
//                 const int p = left[current];
//                 const int q = right[current];
//
//                 expandEndpoint(p, distance);
//                 expandEndpoint(q, distance);
//
//                 /*
//                  * Report every undiscovered link f satisfying
//                  *
//                  *     p < left[f] < q < right[f].
//                  */
//                 const int firstLeftPosition = upperLeft(p);
//                 const int lastLeftPosition = lowerLeft(q);
//
//                 while (true)
//                 {
//                     const int position = maxRightReporter.findPosition(firstLeftPosition, lastLeftPosition, q);
//
//                     if (position == -1)
//                     {
//                         break;
//                     }
//
//                     discover(byLeft[position], distance);
//                 }
//
//                 /*
//                  * Report every undiscovered link f satisfying
//                  *
//                  *     left[f] < p < right[f] < q.
//                  */
//                 const int firstRightPosition = upperRight(p);
//                 const int lastRightPosition = lowerRight(q);
//
//                 while (true)
//                 {
//                     const int position = minLeftReporter.findPosition(firstRightPosition, lastRightPosition, p);
//
//                     if (position == -1)
//                     {
//                         break;
//                     }
//
//                     discover(byRight[position], distance);
//                 }
//             }
//
//             for (int position = groupBegin; position < groupEnd; ++position)
//             {
//                 const int id = testOrder[position];
//                 const int target = targetOf(id, source);
//
//                 if (endpointDistance[target] < weight[id])
//                 {
//                     removable[id] = 1;
//                 }
//             }
//         }
//
//         groupBegin = groupEnd;
//     }
//
//     std::vector<int> answer;
//
//     for (int id = 0; id < numberOfLinks; ++id)
//     {
//         if (removable[id])
//         {
//             answer.push_back(id);
//         }
//     }
//
//     return answer;
// }
// //