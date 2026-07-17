#pragma once

#include <string>
#include <array>
#include <string_view>
#include <stdexcept>
#include <unordered_map>
#include <functional>
#include <memory>
#include <vector>

#include "algorithm_runner.hpp"

using ParamMap = std::unordered_map<std::string, std::string>;

using AlgorithmFactory = std::function<std::unique_ptr<AlgorithmRunner>(const ParamMap&)>;

struct AlgorithmEntry
{
    std::string name;
    std::string description;
    AlgorithmFactory factory;
};

class AlgorithmRegistry
{
public:
    void add(AlgorithmEntry entry);
    void add(std::string name, std::string description, AlgorithmFactory factory);

    std::unique_ptr<AlgorithmRunner> create(const std::string_view name, const ParamMap& params) const;

    std::vector<std::string> names() const;

private:
    std::unordered_map<std::string, AlgorithmEntry> m_entries;
};

AlgorithmRegistry& get_global_registry();
inline AlgorithmRegistry& global_registry = get_global_registry();
