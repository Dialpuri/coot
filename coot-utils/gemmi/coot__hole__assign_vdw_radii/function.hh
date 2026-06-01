#pragma once

#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include <vector>
#include <string>
#include <map>
#include <utility>

namespace coot {
namespace hole {

// Pad a gemmi atom name (trimmed) to 4-char PDB format
inline std::string pad_atom_name(const std::string& name) {
    if (name.size() >= 4) return name.substr(0, 4);
    std::string padded = name;
    // First char is space if name starts with a letter (for names like "N" -> " N  ")
    // Standard PDB format: element-containing atoms have space-prefixed names
    // e.g. "N" -> " N  ", "CA" -> " CA ", "CB" -> " CB ", "C" -> " C  "
    // Multi-char atom names like "HG2" -> "HG2 " are not padded on the left
    // Heuristic: if name has only 1 char and it's a letter, pad on left
    if (name.size() == 1) {
        padded = " " + name + "  ";
    } else {
        while (padded.size() < 4) padded += ' ';
    }
    return padded;
}

struct vdw_radii_result {
    // chain_idx -> res_idx -> atom_idx -> radius
    std::vector<std::vector<std::vector<double>>> radii;

    size_t count() const {
        size_t n = 0;
        for (const auto& chain_radii : radii)
            for (const auto& res_radii : chain_radii)
                n += res_radii.size();
        return n;
    }

    double sum() const {
        double s = 0;
        for (const auto& chain_radii : radii)
            for (const auto& res_radii : chain_radii)
                for (double r : res_radii)
                    s += r;
        return s;
    }
};

vdw_radii_result assign_vdw_radii_gemmi(
    gemmi::Model& model,
    const coot::protein_geometry& geom);

} // namespace hole
} // namespace coot