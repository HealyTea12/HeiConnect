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

struct AlgorithmParameter
{
    std::string name;
    std::string default_value;
    std::string description;
};

struct AlgorithmEntry
{
    std::string name;
    std::string description;
    std::vector<AlgorithmParameter> parameters;
    AlgorithmFactory factory;
};

class AlgorithmRegistry
{
public:
    void add(AlgorithmEntry entry);
    void add(std::string name, std::string description, AlgorithmFactory factory);
    void add(
        std::string name,
        std::string description,
        std::vector<AlgorithmParameter> parameters,
        AlgorithmFactory factory);

    std::unique_ptr<AlgorithmRunner> create(const std::string_view name, const ParamMap& params) const;
    const AlgorithmEntry* find(std::string_view name) const;

    std::vector<std::string> names() const;

private:
    std::unordered_map<std::string, AlgorithmEntry> m_entries;
};

AlgorithmRegistry& get_global_registry();
inline AlgorithmRegistry& global_registry = get_global_registry();
