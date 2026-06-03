#pragma once

#include "geometry/protein-geometry.hh"
#include <memory>
#include <string>
#include <vector>
#include <map>

// Replacement for mmdb::math::Graph, mmdb::math::Vertex, mmdb::math::Edge
struct dict_graph_vertex_t {
    std::string element;
    std::string name;
};

struct dict_graph_edge_t {
    unsigned int v1;  // 1-indexed
    unsigned int v2;  // 1-indexed
    int bond_type;
};

struct dict_graph_t {
    std::vector<dict_graph_vertex_t> vertices;
    std::vector<dict_graph_edge_t> edges;

    unsigned int GetNofVertices() const { return static_cast<unsigned int>(vertices.size()); }
    unsigned int GetNofEdges() const    { return static_cast<unsigned int>(edges.size()); }
};

namespace coot {

// Free-function gemmi port of dictionary_residue_restraints_t::make_graph
std::unique_ptr<dict_graph_t> make_graph_gemmi(const dictionary_residue_restraints_t &drr, bool use_hydrogens);

} // namespace coot