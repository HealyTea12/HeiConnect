#pragma once
#include <chrono>
#include <string>
#include <iostream>
#include <filesystem>
#include <optional>
#include <fstream>

class Timer
{
public:
    Timer(const std::string &name, const std::optional<std::filesystem::path> &output_path = std::nullopt)
        : m_name(name), m_start(std::chrono::high_resolution_clock::now()), m_output_path(output_path)
    {
    }

    ~Timer()
    {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - m_start).count();
        std::cout << "Timer [" << m_name << "]: " << duration << " ms." << std::endl;
        if (m_output_path.has_value())
        {
            auto output_path = m_output_path.value();
            std::ofstream file{output_path.string()};
            if (file.is_open())
            {
                file << m_name << ": " << duration << " ms." << std::endl;
            }
        }
    }

private:
    std::string m_name;
    std::chrono::high_resolution_clock::time_point m_start;
    std::optional<std::filesystem::path> m_output_path;
};