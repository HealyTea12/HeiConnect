#include <filesystem>
#include <fstream>
#include <sstream>

#include "gtest/gtest.h"

#include "HeiConnect/visualization/dot_writer.hpp"

TEST(DotWriterTest, WritesOriginalAndLinkGraphsWithDifferentColors)
{
    const WeightedCRFGraph<> original_graph{{0, 1, 3, 4}, {1, 0, 2, 1}, {1.0, 1.0, 1.0, 1.0}};
    const WeightedCRFGraph<> link_graph{{0, 1, 1, 2}, {2, 0}, {2.0, 2.0}};
    const auto output_path = std::filesystem::temp_directory_path() / "heiconnect_dot_writer_test.dot";

    HeiConnect::visualization::write_to_dot(original_graph, link_graph, output_path);

    std::ifstream output{output_path};
    std::stringstream contents;
    contents << output.rdbuf();

    EXPECT_EQ(
        contents.str(),
        "graph G {\n"
        "    0;\n"
        "    1;\n"
        "    2;\n"
        "    0 -- 1 [color=black];\n"
        "    1 -- 2 [color=black];\n"
        "    0 -- 2 [color=pink];\n"
        "}\n");

    std::filesystem::remove(output_path);
}

TEST(DotWriterTest, RejectsGraphsWithDifferentNumbersOfVertices)
{
    const WeightedCRFGraph<> original_graph{{0, 0}, {}, {}};
    const WeightedCRFGraph<> link_graph{{0, 0, 0}, {}, {}};

    EXPECT_THROW(
        HeiConnect::visualization::write_to_dot(original_graph, link_graph, "unused.dot"),
        std::invalid_argument);
}
