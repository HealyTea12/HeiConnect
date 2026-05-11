#pragma once
#include <chrono>
#include <string>
#include <iostream>
#include <filesystem>
#include <optional>
#include <unordered_map>
#include <fstream>
#include <omp.h>

class Timer
{
public:
    Timer(const std::string &name, const std::optional<std::filesystem::path> &output_path = std::nullopt)
        : m_name(name), m_start(std::chrono::high_resolution_clock::now()), m_output_path(output_path)
    {
    }

    void add_checkpoint(std::string name)
    {
        auto now = std::chrono::high_resolution_clock::now();
        m_checkpoints.push_back({name, now});
    }

    ~Timer()
    {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - m_start).count();
        if (m_output_path.has_value())
        {
            auto output_path = m_output_path.value();
            std::ofstream file{output_path.string(), std::ios::app};
            if (file.is_open())
            {
                file << m_name << ": " << duration << " ms." << std::endl;
                if (!m_checkpoints.empty())
                {
                    auto previous = m_start;
                    file << "{\n";
                    for (const auto &checkpoint : m_checkpoints)
                    {
                        auto checkpoint_duration = std::chrono::duration_cast<std::chrono::milliseconds>(checkpoint.second - previous).count();
                        file << " " << checkpoint.first << ": " << checkpoint_duration << " ms." << std::endl;
                        previous = checkpoint.second;
                    }
                    file << "}\n";
                }
            }
        }
    }

private:
    std::string m_name;
    std::chrono::high_resolution_clock::time_point m_start;
    std::optional<std::filesystem::path> m_output_path;
    std::vector<std::pair<std::string, std::chrono::time_point<std::chrono::high_resolution_clock>>> m_checkpoints{};
};

inline double get_time_s()
{
    return omp_get_wtime();
}