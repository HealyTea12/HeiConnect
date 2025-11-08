#pragma once

#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include "set_cover/graph.hpp"

WeightedCRFGraph read_from_file(const std::filesystem::path &path);

WeightedCRFGraph read_from_file_graphML(const std::filesystem::path &path);