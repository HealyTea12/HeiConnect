#include <gtest/gtest.h>

#include "HeiConnect/conn_aug/reducers/common.hpp"
#include "HeiConnect/data_structures/graph_utils.hpp"

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
