#pragma once

#include <sys/wait.h>
#include <sys/resource.h>
#include <iostream>
#include <fstream>
#include <chrono>

#include "HeiConnect/tools/timer.hpp"

// Helper to run a function in an isolated process and measure its peak memory usage in pages
// linux specific
template <class ChildFn>
long run_isolated_and_measure_memory_usage(ChildFn child_function)
{
    pid_t pid = fork();
    if (pid == -1)
    {
        throw std::runtime_error("Failed to fork process...");
    }
    else if (pid == 0)
    {
        try
        {
            child_function();
            std::exit(0);
        }
        catch (const std::exception &e)
        {
            std::exit(1);
        }
    }
    else
    {
        int status;
        struct rusage rusage{};
        pid_t wpid = wait4(pid, &status, 0, &rusage);

        const long page_size = sysconf(_SC_PAGESIZE);
        const long peak_kb = rusage.ru_maxrss;
        const long peak_bytes = peak_kb * 1024L;
        const long peak_pages = (peak_bytes + page_size - 1) / page_size;
        return peak_pages;
    }
}

// Helper to write to both stdout and output file
inline void log_to_file_and_stdout(const std::string &message, const std::filesystem::path &output_file)
{
    std::cout << message << std::endl;
    std::ofstream ofs{output_file.string(), std::ios::app};
    ofs << message << std::endl;
}

inline void log_separator(const std::filesystem::path &output_file)
{
    std::cout << "----------------------------------------" << std::endl;
    std::ofstream ofs{output_file.string(), std::ios::app};
    ofs << "----------------------------------------" << std::endl;
}
