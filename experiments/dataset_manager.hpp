#include <filesystem>
#include <stdexcept>

class TemporaryDir
{
    TemporaryDir(std::filesystem::path dir) : m_base_dir(dir)
    {
        if (std::filesystem::exists(dir))
            throw std::runtime_error("Directory already exists");
        std::filesystem::create_directory(dir);
    };
    ~TemporaryDir()
    {
        std::filesystem::remove(m_base_dir);
    }

private:
    std::filesystem::path m_base_dir;
};