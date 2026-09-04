#include <gtest/gtest.h>

#include "HeiConnect/conn_aug/reducers/common.hpp"
#include "HeiConnect/conn_aug/reducers/cycle_reducer.hpp"
#include "HeiConnect/data_structures/intersection_index/intersection_tree.hpp"
#include "HeiConnect/data_structures/graph_utils.hpp"

#include <array>
#include <limits>
#include <random>

static std::vector<std::vector<int>> naive_cycle_distance_closure(std::vector<std::vector<int>> distance)
{
    const size_t n = distance.size();
    bool changed = true;
    while (changed)
    {
        changed = false;
        auto relax = [&](size_t u, size_t v, int candidate) {
            if (candidate < distance[u][v])
            {
                distance[u][v] = candidate;
                distance[v][u] = candidate;
                changed = true;
            }
        };

        for (size_t x = 0; x < n; ++x)
        {
            for (size_t y = 0; y < n; ++y)
            {
                for (size_t z = 0; z < n; ++z)
                {
                    relax(x, z, distance[x][y] + distance[y][z]);
                }
            }
        }

        for (size_t a = 0; a < n; ++a)
        {
            for (size_t b = a + 1; b < n; ++b)
            {
                for (size_t c = b + 1; c < n; ++c)
                {
                    for (size_t d = c + 1; d < n; ++d)
                    {
                        const int candidate = distance[a][c] + distance[b][d];
                        relax(a, b, candidate);
                        relax(b, c, candidate);
                        relax(c, d, candidate);
                        relax(a, d, candidate);
                    }
                }
            }
        }
    }
    return distance;
}

TEST(ConnAugReducers, MergesCrossingContractionsOnCycle)
{
    auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr({{0, 2, 1.0}, {1, 3, 1.0}}, 4);
    const std::vector<size_t> selected_links{0, 1};
    UnionFind uf(graph.num_vertices());

    const auto stats = add_links_to_union_find(graph, link_graph, uf, selected_links);

    EXPECT_EQ(uf.find(0), uf.find(1));
    EXPECT_EQ(uf.find(0), uf.find(2));
    EXPECT_EQ(uf.find(0), uf.find(3));
    EXPECT_EQ(stats.total_merged_nodes, 3);
}

TEST(ConnAugReducers, KeepsNonCrossingContractionsSeparateOnCycle)
{
    auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr({{0, 2, 1.0}}, 4);
    const std::vector<size_t> selected_links{0};
    UnionFind uf(graph.num_vertices());

    const auto stats = add_links_to_union_find(graph, link_graph, uf, selected_links);

    EXPECT_EQ(uf.find(0), uf.find(2));
    EXPECT_NE(uf.find(0), uf.find(1));
    EXPECT_NE(uf.find(0), uf.find(3));
    EXPECT_EQ(stats.total_merged_nodes, 1);
}

TEST(ConnAugReducers, KeepsNestedContractionsSeparateOnCycle)
{
    auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr({{0, 3, 1.0}, {1, 2, 1.0}}, 4);
    const std::vector<size_t> selected_links{0, 1};
    UnionFind uf(graph.num_vertices());

    const auto stats = add_links_to_union_find(graph, link_graph, uf, selected_links);

    EXPECT_EQ(uf.find(0), uf.find(3));
    EXPECT_EQ(uf.find(1), uf.find(2));
    EXPECT_NE(uf.find(0), uf.find(1));
    EXPECT_EQ(stats.total_merged_nodes, 2);
}

TEST(ConnAugReducers, FrozenStackMergesCrossingContractionsOnCycle)
{
    auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr({{0, 2, 1.0}, {1, 3, 1.0}}, 4);
    const std::vector<size_t> selected_links{0, 1};
    UnionFind uf(graph.num_vertices());

    const auto stats = add_links_to_union_find_frozen_stack(graph, link_graph, uf, selected_links);

    EXPECT_EQ(uf.find(0), uf.find(1));
    EXPECT_EQ(uf.find(0), uf.find(2));
    EXPECT_EQ(uf.find(0), uf.find(3));
    EXPECT_EQ(stats.total_merged_nodes, 3);
}

TEST(ConnAugReducers, FrozenStackKeepsNonCrossingContractionsSeparateOnCycle)
{
    auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr({{0, 2, 1.0}}, 4);
    const std::vector<size_t> selected_links{0};
    UnionFind uf(graph.num_vertices());

    const auto stats = add_links_to_union_find_frozen_stack(graph, link_graph, uf, selected_links);

    EXPECT_EQ(uf.find(0), uf.find(2));
    EXPECT_NE(uf.find(0), uf.find(1));
    EXPECT_NE(uf.find(0), uf.find(3));
    EXPECT_EQ(stats.total_merged_nodes, 1);
}

TEST(ConnAugReducers, FrozenStackKeepsNestedContractionsSeparateOnCycle)
{
    auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr({{0, 3, 1.0}, {1, 2, 1.0}}, 4);
    const std::vector<size_t> selected_links{0, 1};
    UnionFind uf(graph.num_vertices());

    const auto stats = add_links_to_union_find_frozen_stack(graph, link_graph, uf, selected_links);

    EXPECT_EQ(uf.find(0), uf.find(3));
    EXPECT_EQ(uf.find(1), uf.find(2));
    EXPECT_NE(uf.find(0), uf.find(1));
    EXPECT_EQ(stats.total_merged_nodes, 2);
}

TEST(ConnAugReducers, MaterializedContractionsPreserveParallelLinks)
{
    const auto graph = WeightedCRFGraph<>::vec_links_to_csr(
        {{0, 1, 1.0}, {1, 0, 1.0}, {1, 2, 1.0}, {2, 1, 1.0}},
        3);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr(
        {{0, 2, 2.0}, {0, 1, 1.0}, {1, 2, 3.0}},
        3);
    UnionFind uf(graph.num_vertices());
    uf.unite(0, 1);

    auto [contracted_graph, contracted_link_graph, original_link_ids] =
        materialize_contractions_preserving_links(graph, link_graph, uf);

    EXPECT_EQ(contracted_graph.num_vertices(), 2);
    ASSERT_EQ(contracted_link_graph.num_edges(), 2);
    EXPECT_EQ(contracted_link_graph.graph.edges[0], contracted_link_graph.graph.edges[1]);
    EXPECT_EQ(original_link_ids, (std::vector<size_t>{0, 2}));
}

TEST(ConnAugReducers, IntersectionTreeMatchesBaselineCycleReduction)
{
    std::vector<std::tuple<int, int, int>> links;
    constexpr int cycle_size = 12;
    for (int u = 0; u < cycle_size; ++u)
    {
        for (int v = u + 1; v < cycle_size; ++v)
        {
            links.emplace_back(u, v, (u + v) % cycle_size + 1);
        }
    }

    const BaselineIntersectionIdx<0> baseline;
    const IntersectionTreeIdx<0> intersection_tree;
    const WeightedIntersectionTreeIdx<0> weighted_intersection_tree;
    const auto baseline_result = cycle_domination_baseline(links, cycle_size, baseline);
    EXPECT_EQ(
        baseline_result,
        cycle_domination_baseline(links, cycle_size, intersection_tree));
    EXPECT_EQ(
        baseline_result,
        cycle_domination_baseline(links, cycle_size, weighted_intersection_tree));
}

TEST(ConnAugReducers, DistinguishesTouchingCrossingAndContainedIntervals)
{
    EXPECT_EQ(interval_overlap({0, 2}, {3, 4}), -1);
    EXPECT_EQ(interval_overlap({0, 2}, {2, 4}), 0);
    EXPECT_EQ(interval_overlap({0, 2}, {1, 3}), 1);
    EXPECT_EQ(interval_overlap({0, 4}, {1, 3}), 2);

    EXPECT_TRUE(intervals_touch_or_cross({0, 2}, {2, 4}));
    EXPECT_TRUE(intervals_touch_or_cross({0, 2}, {1, 3}));
    EXPECT_TRUE(intervals_touch_or_cross({0, 4}, {0, 3}));
    EXPECT_FALSE(intervals_touch_or_cross({0, 2}, {3, 4}));
    EXPECT_FALSE(intervals_touch_or_cross({0, 4}, {1, 3}));
}

TEST(ConnAugReducers, IncrementalBaselineStartsEmptyAndStopsAtExclusiveLevel)
{
    const std::vector<IntersectionRecord> possible_records{
        {{0, 2}, 0},
        {{1, 3}, 0},
        {{2, 4}, 0}};
    const BaselineIntersectionIdx<2> index_type;
    auto index = index_type.makeEmpty(possible_records);
    std::vector<size_t> intersections;

    index->forEachIntersection(
        [&](size_t id, IntersectionInterval) { intersections.push_back(id); },
        {0, 4},
        3);
    EXPECT_TRUE(intersections.empty());

    EXPECT_EQ(index->addInterval({{0, 2}, 0}), 0);
    EXPECT_EQ(index->addInterval({{1, 3}, 1}), 1);
    EXPECT_EQ(index->addInterval({{2, 4}, 2}), 2);
    index->forEachIntersection(
        [&](size_t id, IntersectionInterval) { intersections.push_back(id); },
        {0, 4},
        2);

    EXPECT_EQ(intersections, (std::vector<size_t>{0}));
    EXPECT_EQ(index->emit_metrics().candidates_inspected, 2);
}

TEST(ConnAugReducers, DistanceClosureAppliesTriangleRule)
{
    const std::vector<std::vector<int>> distance{{0, 2, 10}, {2, 0, 3}, {10, 3, 0}};
    const auto closed = cycle_distance_closure(distance);

    EXPECT_EQ(closed[0][2], 5);
}

TEST(ConnAugReducers, DistanceClosureAppliesCrossingRule)
{
    const std::vector<std::vector<int>> distance{
        {0, 10, 2, 10},
        {10, 0, 10, 3},
        {2, 10, 0, 10},
        {10, 3, 10, 0}};
    const auto closed = cycle_distance_closure(distance);

    EXPECT_EQ(closed[0][1], 5);
    EXPECT_EQ(closed[1][2], 5);
    EXPECT_EQ(closed[2][3], 5);
    EXPECT_EQ(closed[0][3], 5);
}

TEST(ConnAugReducers, DistanceClosureRecordsIndexIterationsAtMetricsLevelTwo)
{
    const std::vector<std::vector<int>> distance{
        {0, 10, 2, 10},
        {10, 0, 10, 3},
        {2, 10, 0, 10},
        {10, 3, 10, 0}};
    const BaselineIntersectionIdx<2> intersection_index;
    CycleReductionMetrics metrics;

    cycle_distance_closure<2>(distance, intersection_index, &metrics);

    constexpr size_t number_of_pairs = 6;
    constexpr size_t worst_case_iterations = number_of_pairs * (number_of_pairs + 1) / 2;
    EXPECT_GT(metrics.intersection_index.candidates_inspected, 0);
    EXPECT_LE(metrics.intersection_index.candidates_inspected, worst_case_iterations);
}

TEST(ConnAugReducers, DistanceClosureTightensMaximumAfterRelaxation)
{
    const std::vector<std::vector<int>> distance{{0, 1, 5}, {1, 0, 1}, {5, 1, 0}};
    const BaselineIntersectionIdx<2> intersection_index;
    CycleReductionMetrics metrics;

    const auto closed = cycle_distance_closure<2>(distance, intersection_index, &metrics);

    EXPECT_EQ(closed[0][2], 2);
    EXPECT_EQ(metrics.priority_queue_pops, 2);
    EXPECT_EQ(metrics.intersection_index.candidates_inspected, 1);
    EXPECT_EQ(metrics.termination_by_cutoff, 1);
}

TEST(ConnAugReducers, DistanceClosureDoesNotApplyCrossingRuleToContainment)
{
    const std::vector<std::vector<int>> distance{
        {0, 100, 100, 1},
        {100, 0, 1, 100},
        {100, 1, 0, 100},
        {1, 100, 100, 0}};
    const auto closed = cycle_distance_closure(distance);

    EXPECT_EQ(closed, distance);
}

TEST(ConnAugReducers, DistanceClosureMatchesNaiveFixedPoint)
{
    std::mt19937 random_engine(246810);
    std::uniform_int_distribution<int> link_weight(1, 100);
    for (size_t n = 2; n <= 8; ++n)
    {
        for (int instance = 0; instance < 100; ++instance)
        {
            std::vector<std::vector<int>> distance(n, std::vector<int>(n));
            for (size_t u = 0; u < n; ++u)
            {
                for (size_t v = u + 1; v < n; ++v)
                {
                    distance[u][v] = link_weight(random_engine);
                    distance[v][u] = distance[u][v];
                }
            }

            const auto expected = naive_cycle_distance_closure(distance);
            EXPECT_EQ(cycle_distance_closure(distance), expected);
            EXPECT_EQ(cycle_distance_closure<IntersectionTreeIdx<0>>(distance), expected);
            EXPECT_EQ(cycle_distance_closure<WeightedIntersectionTreeIdx<0>>(distance), expected);
        }
    }
}

TEST(ConnAugReducers, GlobalSweepMatchesBaselineOnFullLinkGraph)
{
    std::vector<std::tuple<int, int, int>> links;
    constexpr int cycle_size = 12;
    for (int u = 0; u < cycle_size; ++u)
    {
        for (int v = u + 1; v < cycle_size; ++v)
        {
            links.emplace_back(u, v, (u + v) % cycle_size + 1);
        }
    }

    const auto expected = cycle_domination_baseline(links, cycle_size);
    const IntersectionTreeIdx<0> intersection_tree;
    const WeightedIntersectionTreeIdx<0> weighted_intersection_tree;

    EXPECT_EQ(cycle_domination_global_sweep(links, cycle_size), expected);
    EXPECT_EQ(cycle_domination_global_sweep(links, cycle_size, intersection_tree), expected);
    EXPECT_EQ(cycle_domination_global_sweep(links, cycle_size, weighted_intersection_tree), expected);
}

TEST(ConnAugReducers, GlobalSweepMatchesBaselineOnRandomCompleteTables)
{
    std::mt19937 random_engine(123456);
    std::uniform_int_distribution<int> link_weight(1, 1000);
    const WeightedIntersectionTreeIdx<0> weighted_intersection_tree;
    for (int cycle_size = 3; cycle_size <= 10; ++cycle_size)
    {
        for (int instance = 0; instance < 100; ++instance)
        {
            std::vector<std::tuple<int, int, int>> links;
            for (int u = 0; u < cycle_size; ++u)
            {
                for (int v = u + 1; v < cycle_size; ++v)
                {
                    links.emplace_back(u, v, link_weight(random_engine));
                }
            }

            const auto expected = cycle_domination_baseline(links, cycle_size, weighted_intersection_tree);
            EXPECT_EQ(cycle_domination_global_sweep(links, cycle_size), expected);
            EXPECT_EQ(cycle_domination_global_sweep(links, cycle_size, weighted_intersection_tree), expected);
        }
    }
}

TEST(ConnAugReducers, GlobalSweepMatchesBaselineExhaustivelyOnFourVertices)
{
    const std::array<std::pair<int, int>, 6> endpoints{
        std::pair{0, 1},
        std::pair{0, 2},
        std::pair{0, 3},
        std::pair{1, 2},
        std::pair{1, 3},
        std::pair{2, 3}};

    constexpr int assignments = 4 * 4 * 4 * 4 * 4 * 4;
    for (int assignment = 0; assignment < assignments; ++assignment)
    {
        int encoded = assignment;
        std::vector<std::tuple<int, int, int>> links;
        for (const auto [u, v] : endpoints)
        {
            const int weight = encoded % 4 + 1;
            encoded /= 4;
            links.emplace_back(u, v, weight);
        }

        EXPECT_EQ(cycle_domination_global_sweep(links, 4), cycle_domination_baseline(links, 4));
    }
}

TEST(ConnAugReducers, GlobalSweepInitializesEveryPairState)
{
    const std::vector<std::tuple<int, int, int>> links{
        {0, 1, 1}, {0, 2, 8}, {0, 3, 9}, {1, 2, 2}, {1, 3, 7}, {2, 3, 3}};
    CycleReductionMetrics metrics;

    cycle_domination_global_sweep<2>(links, 4, &metrics);

    EXPECT_EQ(metrics.sources, 6);
}

TEST(ConnAugReducers, GlobalSweepKeepsStrictEquality)
{
    const std::vector<std::tuple<int, int, int>> links{
        {0, 1, 4}, {0, 2, 2}, {0, 3, 100}, {1, 2, 100}, {1, 3, 2}, {2, 3, 100}};
    const WeightedIntersectionTreeIdx<0> weighted_intersection_tree;
    const auto expected = cycle_domination_baseline(links, 4, weighted_intersection_tree);
    const auto actual = cycle_domination_global_sweep(links, 4, weighted_intersection_tree);

    EXPECT_EQ(actual, expected);
    EXPECT_EQ(std::find(actual.begin(), actual.end(), 0), actual.end());
}

TEST(ConnAugReducers, GlobalSweepHandlesWeightCutoffWithoutOverflow)
{
    const auto maximum_weight = std::numeric_limits<uint64_t>::max();
    const std::vector<std::tuple<int, int, uint64_t>> links{
        {0, 1, maximum_weight},
        {0, 2, maximum_weight - 2},
        {0, 3, maximum_weight},
        {1, 2, maximum_weight},
        {1, 3, 2},
        {2, 3, maximum_weight}};
    const auto expected = cycle_domination_baseline(links, 4);
    const auto actual = cycle_domination_global_sweep(links, 4);

    EXPECT_EQ(actual, expected);
    EXPECT_EQ(std::find(actual.begin(), actual.end(), 0), actual.end());
}

TEST(ConnAugReducers, GlobalSweepUsesTheMaximumDistanceCutoff)
{
    const std::vector<std::tuple<int, int, int>> links{{0, 1, 1}, {1, 2, 2}, {0, 2, 10}};
    CycleReductionMetrics metrics;

    cycle_domination_global_sweep<2>(links, 3, &metrics);

    EXPECT_EQ(metrics.sources, 3);
    EXPECT_EQ(metrics.termination_by_completion, 0);
    EXPECT_EQ(metrics.termination_by_cutoff + metrics.termination_by_empty_queue, 1);
}

TEST(ConnAugReducers, RecordsDetailedCycleReductionMetrics)
{
    const std::vector<std::tuple<int, int, int>> links{{0, 2, 1}, {1, 3, 2}, {0, 3, 3}};
    const WeightedIntersectionTreeIdx<2> intersection_tree;
    CycleReductionMetrics metrics;

    cycle_domination_baseline<2>(links, 4, intersection_tree, &metrics);

    EXPECT_EQ(metrics.sources, 4);
    EXPECT_GT(metrics.priority_queue_pops, 0);
    EXPECT_LE(metrics.priority_queue_pops, metrics.possible_priority_queue_pops);
    EXPECT_GT(metrics.intersection_index.queries, 0);
    EXPECT_GT(
        metrics.intersection_index.candidates_inspected + metrics.intersection_index.subtrees_pruned_by_level,
        0);
    EXPECT_GT(metrics.maximum_pops_per_link, 0);
    EXPECT_EQ(
        metrics.intersection_index.callbacks,
        metrics.intersection_candidates_enqueued + metrics.intersection_candidates_already_explored +
            metrics.intersection_candidates_rejected_by_cutoff);
    EXPECT_EQ(metrics.intersection_candidates_already_explored, 0);
    EXPECT_EQ(metrics.intersection_candidates_rejected_by_cutoff, 0);
}
