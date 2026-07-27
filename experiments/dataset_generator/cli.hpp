#pragma once

#include "config.hpp"

#include <optional>

namespace dataset_generator
{

std::optional<Config> parse_config(int argc, char **argv);

} // namespace dataset_generator
