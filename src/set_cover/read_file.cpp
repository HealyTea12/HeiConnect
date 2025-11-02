#include "set_cover/read_file.hpp"

// construct CRF graph from xml file
WeightedCRFGraph read_from_file(const std::filesystem::path &path)
{
    size_t n, m;
    if (!std::filesystem::exists(path))
    {
        throw std::runtime_error("File not found: " + path.string());
    }
    std::ifstream file{path};
    std::string line;
    do
    {
        std::getline(file, line);
    } while (line[0] == '%');
    std::istringstream iss(line);
    iss >> n >> m;
    int weight_fmt = 0, n_vertex_weights = 0, tmp;
    if (!iss.eof())
    {
        iss >> weight_fmt;
    }
    if (weight_fmt >= 10)
    {
        if (iss.eof())
        {
            n_vertex_weights = 1;
        }
        else
        {
            iss >> n_vertex_weights;
        }
    }
    bool has_edge_weights = weight_fmt % 2 == 1;
    unsigned int u = 0, v = 0;
    std::vector<size_t> vertices;
    std::vector<size_t> edges;
    std::vector<double> weights;
    while (std::getline(file, line))
    {
        if (line[0] == '%')
        {
            continue;
        }
        std::istringstream iss(line);
        for (int i = 0; i < n_vertex_weights; ++i)
        {
            iss >> tmp;
        }
        while ((iss >> v))
        {
            v--; // zero-based indexing
            if (u < v)
            {
                edges.push_back(v);
                // add_edge(u, v);
                if (has_edge_weights)
                {
                    iss >> tmp;
                    weights.push_back(static_cast<double>(tmp));
                }
            }
        }
        vertices.push_back(edges.size());
        ++u;
    }
    return WeightedCRFGraph{{vertices, edges}, weights};
}