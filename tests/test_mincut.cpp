#include <gtest/gtest.h>

#include "HeiConnect/min_cut/simple_mincut.hpp"

TEST(MinCutTest, SimpleGraph)
{
    // Create a simple weighted graph
    WeightedCRFGraph graph(
        CRFGraph(
            {0, 2, 4, 6},      // vertices
            {1, 2, 0, 2, 0, 1} // edges
            ),
        {1.0, 2.0, 1.0, 3.0, 2.0, 3.0} // weights
    );

    double min_cut = global_mincut_simple(graph);
    EXPECT_DOUBLE_EQ(min_cut, 3.0);
}