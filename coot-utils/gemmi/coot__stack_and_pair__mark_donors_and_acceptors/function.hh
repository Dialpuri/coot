#pragma once

#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include <map>
#include <string>
#include <vector>

namespace coot {
namespace stack_and_pair {

// Pad a gemmi atom name to PDB-format 4 characters (right-justified,
// padded on the LEFT with spaces), matching what MMDB and protein_geometry
// dictionary keys expect.
static inline std::string pad_atom_name(const std::string& name) {
    if (name.size() == 1) return " " + name + "  ";
    if (name.size() == 2) return " " + name + " ";
    if (name.size() == 3) return " " + name;
    return name;  // already 4+ chars
}

// gemmi port of mark_donors_and_acceptors:
// For each atom (given as CRA), look up or compute the H-bond type
// using protein_geometry::get_h_bond_type, caching results by
// atom_name+res_name. Returns a parallel vector of hb_type ints.
//
// MMDB original stored results in per-atom UDD and returned the
// UDD handle.  gemmi has no UDD, so the vector replaces that
// mechanism.
std::vector<int> mark_donors_and_acceptors_gemmi(
    const std::vector<gemmi::CRA>& atoms,
    const protein_geometry& geom,
    int imol_enc)
{
    std::vector<int> hb_types;
    std::map<std::string, int> cache;

    for (const auto& cra : atoms) {
        if (!cra.atom) continue;

        std::string name = pad_atom_name(cra.atom->name);
        std::string res_name = cra.residue ? cra.residue->name : "";

        std::string key = name + "+" + res_name;
        auto it = cache.find(key);
        int hb_type;
        if (it != cache.end()) {
            hb_type = it->second;
        } else {
            hb_type = geom.get_h_bond_type(name, res_name, imol_enc);
            cache[key] = hb_type;
        }
        hb_types.push_back(hb_type);
    }

    return hb_types;
}

} // namespace stack_and_pair
} // namespace coot