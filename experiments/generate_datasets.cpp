#include "dataset_generator/cli.hpp"
#include "dataset_generator/generators.hpp"

#include <iostream>

int main(int argc, char **argv)
{
    try
    {
        auto config = dataset_generator::parse_config(argc, argv);
        if (!config)
        {
            return 0;
        }

        dataset_generator::generate(*config);
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
