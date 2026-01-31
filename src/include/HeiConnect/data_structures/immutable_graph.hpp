#pragma once
#include <vector>
#include <ranges>
#include <unordered_map>
#include <unordered_set>
#include <concepts>
#include <filesystem>
#include <fstream>

#include "pugixml.hpp"

template <typename NodeID = size_t, typename EdgeID = size_t>
struct CRFGraph
{
    std::vector<NodeID> vertices;
    std::vector<EdgeID> edges;

    CRFGraph(const std::vector<NodeID> &vertices,
             const std::vector<EdgeID> &edges)
        : vertices(vertices), edges(edges) {}
};

template <typename NodeID = size_t, typename EdgeID = size_t, typename WeightType = double>
struct WeightedCRFGraph
{
    CRFGraph<NodeID, EdgeID> graph;
    std::vector<WeightType> weights;

public:
    WeightedCRFGraph(const CRFGraph<NodeID, EdgeID> &graph,
                     const std::vector<WeightType> &weights)
        : graph(graph), weights(weights) {}

    WeightedCRFGraph(std::vector<NodeID> vertices,
                     std::vector<EdgeID> edges,
                     std::vector<WeightType> weights)
        : graph{std::move(vertices), std::move(edges)}, weights{std::move(weights)} {}

    bool is_edge(NodeID u, NodeID v) const
    {
        for (EdgeID e{graph.vertices[u]}; e < graph.vertices[u + 1]; ++e)
        {
            if (graph.edges[e] == v)
                return true;
        }
        return false;
    }

    // Metrics
    size_t num_edges() const
    {
        return graph.edges.size();
    }

    size_t num_vertices() const
    {
        return graph.vertices.size() - 1;
    }

    size_t degree(NodeID u) const
    {
        return graph.vertices[u + 1] - graph.vertices[u];
    }

    double average_degree() const
    {
        return (double)num_edges() / (double)num_vertices();
    }

    size_t max_degree() const
    {
        size_t max_deg = 0;
        for (NodeID u{0}; u < num_vertices(); ++u)
        {
            size_t deg = degree(u);
            if (deg > max_deg)
            {
                max_deg = deg;
            }
        }
        return max_deg;
    }

    size_t min_degree() const
    {
        size_t min_deg = SIZE_MAX;
        for (NodeID u{0}; u < num_vertices(); ++u)
        {
            size_t deg = degree(u);
            if (deg < min_deg)
            {
                min_deg = deg;
            }
        }
        return min_deg;
    }

    // ----

    WeightedCRFGraph make_bidirectional() const
    {
        auto new_vertices = std::vector<size_t>(graph.vertices.size(), 0);
        std::vector<size_t> new_edges{};
        std::vector<double> new_weights{};
        for (size_t u = 0; u < graph.vertices.size() - 1; ++u)
        {
            // add original edges
            for (size_t e = graph.vertices[u]; e < graph.vertices[u + 1]; ++e)
            {
                auto v = graph.edges[e];
                new_edges.emplace_back(v);
                new_weights.emplace_back(weights[e]);
            }
            // add reverse edges
            for (size_t v = 0; v < graph.vertices.size() - 1; ++v)
            {
                for (size_t e = graph.vertices[v]; e < graph.vertices[v + 1]; ++e)
                {
                    if (graph.edges[e] == u)
                    {
                        new_edges.emplace_back(v);
                        new_weights.emplace_back(weights[e]);
                    }
                }
            }
            new_vertices[u + 1] = new_edges.size();
        }
        return WeightedCRFGraph{
            {new_vertices, new_edges},
            new_weights};
    }

    static WeightedCRFGraph read_from_file(const std::filesystem::path &path)
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
        } while (line.length() == 0 || line[0] == '%');
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
        auto vertices = std::vector<size_t>(n + 1, 0);
        auto edges = std::vector<size_t>();
        std::vector<double> weights{};
        while (std::getline(file, line))
        {
            if (line[0] == '%')
            {
                continue;
            }
            std::istringstream iss(line);
            // skip vertex weights
            for (int i = 0; i < n_vertex_weights; ++i)
            {
                iss >> tmp;
            }
            while ((iss >> v))
            {
                v--; // zero-based indexing
                edges.push_back(v);
                if (has_edge_weights)
                {
                    iss >> tmp;
                    weights.push_back(static_cast<double>(tmp));
                }
                else
                {
                    // maybe should warn that unweighted graph
                    weights.push_back(1.0);
                }
            }
            vertices[u + 1] = edges.size();
            ++u;
        }
        return WeightedCRFGraph{{vertices, edges}, weights};
    }

    void write_to_file_metis(const std::filesystem::path &path) const
    {
        std::ofstream file{path};
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file for writing: " + path.string());
        }
        size_t n = graph.vertices.size() - 1;
        size_t m = graph.edges.size() / 2; // undirected graph
        file << n << " " << m << " 1\n";
        for (size_t u = 0; u < n; ++u)
        {
            for (size_t e = graph.vertices[u]; e < graph.vertices[u + 1]; ++e)
            {
                size_t v = graph.edges[e];
                file << v + 1 << " " << static_cast<int>(weights[e]) << " ";
            }
            file << "\n";
        }
    }

    // Construct CSR graph from GraphML file specifically generated by VieCut's mincut cactus output
    // Notes:
    //  - file includes the function mapping the cactus to the original graph
    //  - nodes are zero indexed
    //  - it is an undirected graph, but only one direction is stored
    // TODO: maybe ifdef if we don't care about this file format to get rid of dependency
    static WeightedCRFGraph read_from_file_graphML(const std::filesystem::path &path)
    {
        if (!std::filesystem::exists(path))
            throw std::runtime_error("File not found: " + path.string());
        pugi::xml_document doc;
        pugi::xml_parse_result result = doc.load_file(path.c_str());
        if (!result)
        {
            throw std::runtime_error("Could not open file " + path.string());
        }
        pugi::xml_node xml_graph = doc.child("graphml").child("graph");
        size_t n_nodes = 0;
        size_t n_edges = 0;
        for (auto &node : xml_graph.children("node"))
        {
            ++n_nodes;
        }
        for (auto &edge : xml_graph.children("edge"))
        {
            ++n_edges;
        }
        // for now just use size_t for indices, but need to improve
        auto vertices = std::vector<size_t>(n_nodes + 1, 0);
        auto edges = std::vector<size_t>();
        auto weights = std::vector<double>();
        // parse contained vertices
        // for (auto &node : xml_graph.children("node"))
        //{
        //    size_t node1 = static_cast<size_t>(node.attribute("id").as_ullong()) - 1;
        //    std::string nodes_string =
        //        node.find_child_by_attribute("key", "containedVertices")
        //            .text()
        //            .as_string();

        //    std::stringstream ss(nodes_string);
        //    std::string item;
        //    while (getline(ss, item, ','))
        //    {
        //        size_t node2 = std::stoul(item) - 1;
        //        edges.emplace_back(node2);
        //        graph.map_to_original_graph[cactus_node_id + 1].push_back(nd + 1);
        //        graph.map_to_cactus[nd + 1] = cactus_node_id + 1;
        //    }
        //}

        auto adj_list = std::vector<std::vector<std::pair<size_t, double>>>(n_nodes, std::vector<std::pair<size_t, double>>{});

        // parse edges
        for (auto &edge : xml_graph.children("edge"))
        {
            double weight =
                edge.find_child_by_attribute("key", "weight").text().as_double();
            size_t source = edge.attribute("source").as_ullong();
            size_t target = edge.attribute("target").as_ullong();
            adj_list[source].emplace_back(target, weight);
            adj_list[target].emplace_back(source, weight);
        }
        for (size_t u = 0; u < adj_list.size(); u++)
        {
            for (auto &neighbor : adj_list[u])
            {
                weights.emplace_back(neighbor.second);
                edges.emplace_back(neighbor.first);
            }
            vertices[u + 1] = edges.size();
        }
        return {{vertices, edges}, weights};
    }

    void write_to_file_graphML(const std::filesystem::path &path) const
    {
        pugi::xml_document doc;
        auto decl = doc.append_child(pugi::node_declaration);
        decl.append_attribute("version") = "1.0";
        decl.append_attribute("encoding") = "UTF-8";

        auto graphml = doc.append_child("graphml");
        graphml.append_attribute("xmlns") = "http://graphml.graphdrawing.org/xmlns";
        graphml.append_attribute("xmlns:xsi") = "http://www.w3.org/2001/XMLSchema-instance";
        graphml.append_attribute("xsi:schemaLocation") =
            "http://graphml.graphdrawing.org/xmlns "
            "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd";

        auto graph_node = graphml.append_child("graph");
        graph_node.append_attribute("edgedefault") = "undirected";

        for (size_t u = 0; u < graph.vertices.size() - 1; ++u)
        {
            auto node = graph_node.append_child("node");
            node.append_attribute("id") = std::to_string(u).c_str();
        }

        for (size_t u = 0; u < graph.vertices.size() - 1; ++u)
        {
            for (size_t e = graph.vertices[u]; e < graph.vertices[u + 1]; ++e)
            {
                if (u > graph.edges[e]) // to avoid writing both directions
                    continue;
                auto edge = graph_node.append_child("edge");
                edge.append_attribute("source") = std::to_string(u).c_str();
                edge.append_attribute("target") = std::to_string(graph.edges[e]).c_str();
                auto weight_key = edge.append_child("data");
                weight_key.append_attribute("key") = "weight";
                weight_key.text().set(std::to_string(weights[e]).c_str());
            }
        }
        doc.save_file(path.c_str());
    }

    void write_to_file_graphML(const std::filesystem::path &path,
                               const std::unordered_map<size_t, std::vector<size_t>> &map_to_original_graph) const
    {
        pugi::xml_document doc;
        auto decl = doc.append_child(pugi::node_declaration);
        decl.append_attribute("version") = "1.0";
        decl.append_attribute("encoding") = "UTF-8";

        auto graphml = doc.append_child("graphml");
        graphml.append_attribute("xmlns") = "http://graphml.graphdrawing.org/xmlns";
        graphml.append_attribute("xmlns:xsi") = "http://www.w3.org/2001/XMLSchema-instance";
        graphml.append_attribute("xsi:schemaLocation") =
            "http://graphml.graphdrawing.org/xmlns "
            "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd";

        auto graph_node = graphml.append_child("graph");
        graph_node.append_attribute("edgedefault") = "undirected";

        for (size_t u = 0; u < graph.vertices.size() - 1; ++u)
        {
            auto node = graph_node.append_child("node");
            node.append_attribute("id") = std::to_string(u).c_str();
            auto contained_vertices = node.append_child("data");
            contained_vertices.append_attribute("key") = "containedVertices";
            if (map_to_original_graph.find(u) != map_to_original_graph.end())
            {
                std::string vertices_str;
                for (auto &v : map_to_original_graph.at(u))
                {
                    vertices_str += std::to_string(v) + ",";
                }
                // remove last comma
                if (!vertices_str.empty())
                {
                    vertices_str.pop_back();
                }
                contained_vertices.text().set(vertices_str.c_str());
            }
            else
            {
                contained_vertices.text().set("");
            }
        }

        for (size_t u = 0; u < graph.vertices.size() - 1; ++u)
        {
            for (size_t e = graph.vertices[u]; e < graph.vertices[u + 1]; ++e)
            {
                if (u > graph.edges[e]) // to avoid writing both directions
                    continue;
                auto edge = graph_node.append_child("edge");
                edge.append_attribute("source") = std::to_string(u).c_str();
                edge.append_attribute("target") = std::to_string(graph.edges[e]).c_str();
                auto weight_key = edge.append_child("data");
                weight_key.append_attribute("key") = "weight";
                weight_key.text().set(std::to_string(weights[e]).c_str());
            }
        }
        doc.save_file(path.c_str());
    }

    // Read graph from simple link file format: first line contains n m wf
    // followed by m lines of u v w representing an edge between u and v with weight
    static inline WeightedCRFGraph<> read_from_file_links(const std::filesystem::path &link_file)
    {
        std::ifstream file{link_file};
        size_t n, m, wf;
        file >> n >> m >> wf;
        std::vector<size_t> vertices(n + 1, 0);
        std::vector<size_t> edges;
        std::vector<double> weights;
        size_t u;
        size_t v;
        double w;
        auto adj_list = std::vector<std::vector<std::pair<size_t, double>>>(n);
        for (size_t i = 0; i < m; ++i)
        {
            file >> u >> v >> w;
            --u;
            --v;
            adj_list[u].emplace_back(v, w);
        }
        for (size_t i = 0; i < n; ++i)
        {
            for (const auto &[neighbor, weight] : adj_list[i])
            {
                edges.emplace_back(neighbor);
                weights.emplace_back(weight);
            }
            vertices[i + 1] = edges.size();
        }
        return {{vertices, edges}, weights};
    }

    /// Generate a graph consisting of all non-existing edges in the input graph.
    /// The weights of the new edges are determined by the provided edge_weight function.
    /// @param graph
    /// @param edge_weight
    /// @return
    WeightedCRFGraph generate_links(double (*edge_weight)(size_t u, size_t v)) const
    {
        std::vector<EdgeID> new_vertices{};
        std::vector<NodeID> new_edges{};
        std::vector<double> weights{};
        new_vertices.push_back(0);
        for (size_t u = 0; u < graph.vertices.size() - 1; ++u)
        {
            std::unordered_set<size_t> neighbours{};
            for (size_t idx = graph.vertices[u]; idx < graph.vertices[u + 1]; ++idx)
            {
                neighbours.insert(graph.edges[idx]);
            }
            for (size_t i = u + 1; i < graph.vertices.size() - 1; i++)
            {
                if (i == u)
                    continue;
                if (neighbours.find(i) == neighbours.end())
                {
                    auto weight = edge_weight(u, i);
                    weights.push_back(weight);
                    new_edges.push_back(i);
                }
            }
            new_vertices.push_back(new_edges.size());
        }
        return WeightedCRFGraph{{new_vertices, new_edges}, weights};
    }

    WeightedCRFGraph
    add_links(
        const WeightedCRFGraph &link_graph,
        const std::unordered_set<size_t> &selected_edges) const
    {
        auto new_vertices = std::vector<size_t>(graph.vertices.size(), 0);
        std::vector<size_t> new_edges{};
        std::vector<double> new_weights{};
        for (size_t u = 0; u < graph.vertices.size() - 1; ++u)
        {
            for (size_t e = graph.vertices[u]; e < graph.vertices[u + 1]; ++e)
            {
                auto v = graph.edges[e];
                new_edges.emplace_back(v);
                new_weights.emplace_back(weights[e]);
            }
            for (size_t e = link_graph.graph.vertices[u]; e < link_graph.graph.vertices[u + 1]; ++e)
            {
                if (selected_edges.find(e) != selected_edges.end())
                {
                    auto v = link_graph.graph.edges[e];
                    new_edges.emplace_back(v);
                    new_weights.emplace_back(0.0); // link weights are just the cost of adding them
                }
            }
            new_vertices[u + 1] = new_edges.size();
        }
        return WeightedCRFGraph{
            {new_vertices, new_edges},
            new_weights};
    }

    // will probably get rid of this and construct directly the link
    // graph as a vector of links.
    std::vector<std::tuple<NodeID, NodeID, WeightType>> csr_to_vec_links()
    {
        std::vector<std::tuple<NodeID, NodeID, WeightType>> links = std::vector<std::tuple<NodeID, NodeID, WeightType>>(weights.size());
        for (NodeID u{0}; u < graph.vertices.size() - 1; u++)
        {
            for (EdgeID e{graph.vertices[u]}; e < graph.vertices[u + 1]; e++)
            {
                NodeID v = graph.edges[e];
                WeightType w = weights[e];
                links[e] = {u, v, w};
            }
        }
        return links;
    }

    // TODO: Construction site
    std::vector<size_t> calculate_undirected_indices() const
    {
        std::vector<size_t> uidx = std::vector<size_t>(graph.edges.size(), 0);
        size_t curr_idx{0};
        for (size_t u = 0; u < graph.vertices.size() - 1; ++u)
        {
            for (size_t e = graph.vertices[u]; e < graph.vertices[u + 1]; ++e)
            {
                size_t v = graph.edges[e];
                if (u < v)
                {
                    uidx[e] = curr_idx;
                    // find the reverse edge
                    for (size_t rev_e = graph.vertices[v]; rev_e < graph.vertices[v + 1]; ++rev_e)
                    {
                        if (graph.edges[rev_e] == u)
                        {
                            uidx[rev_e] = curr_idx;
                            break;
                        }
                    }
                    curr_idx++;
                }
            }
        }
        return uidx;
    }

    // TODO: under construction
    // Only works for cactus graphs
    WeightedCRFGraph<> cactus_generate_block_tree(size_t root) const
    {
        auto tin = std::vector<size_t>(graph.vertices.size() - 1, 0);
        auto tout = std::vector<size_t>(graph.vertices.size() - 1, 0);
        auto depth = std::vector<size_t>(graph.vertices.size() - 1, 0);
        auto parent = std::vector<decltype(graph.edges[0])>(graph.vertices.size() - 1, 0);
        size_t timer = 0;
        size_t block_count = 0;
        // dfs
        std::vector<char> visited(graph.vertices.size() - 1, 0);
        std::vector<size_t> stack{root};

        return *this;
    }
};