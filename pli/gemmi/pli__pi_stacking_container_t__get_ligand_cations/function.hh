#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

#include "clipper/core/coords.h"
#include "geometry/protein-geometry.hh"

namespace pli {

namespace {

// Pad atom name to 4 chars in MMDB style (left-padded with spaces)
inline std::string atom_name_to_4c(std::string const &name) {
    std::string trimmed = name;
    size_t a = trimmed.find_first_not_of(' ');
    size_t b = trimmed.find_last_not_of(' ');
    if (a == std::string::npos) return "    ";
    trimmed = trimmed.substr(a, b - a + 1);
    int ilen = static_cast<int>(trimmed.length());
    if (ilen >= 4) return trimmed.substr(0, 4);
    return std::string(4 - ilen, ' ') + trimmed;
}

// Check if an atom is nitrogen: use element info, fall back to name heuristic
inline bool is_nitrogen_atom(gemmi::Atom const &atom) {
    std::string ele = atom.element.name();
    if (!ele.empty()) {
        // Trim trailing whitespace (some CIF readers leave it)
        while (!ele.empty() && ele.back() == ' ') ele.pop_back();
        if (ele == "N") return true;
    }
    // Fallback: if element is unknown (e.g. for non-standard ligands in CIF),
    // check if the atom name starts with 'N' (common convention for N atoms)
    std::string trimmed = atom.name;
    size_t a = trimmed.find_first_not_of(' ');
    if (a != std::string::npos && trimmed[a] == 'N') return true;
    return false;
}

} // anonymous namespace

class pi_stacking_container_t {
public:
    static std::vector<std::pair<std::string, clipper::Coord_orth>>
    get_ligand_cations_gemmi(
        gemmi::Residue const &res_ref,
        coot::dictionary_residue_restraints_t const &monomer_restraints)
    {
        std::vector<std::pair<std::string, clipper::Coord_orth>> v;

        for (gemmi::Atom const &atom : res_ref.atoms) {
            if (!is_nitrogen_atom(atom)) continue;

            // how many bonds does this N have?
            int n_bonds = 0;
            std::string atom_name = atom_name_to_4c(atom.name);

            for (coot::dict_bond_restraint_t const &br : monomer_restraints.bond_restraint) {
                if (br.atom_id_1_4c() == atom_name) {
                    std::string other_atom_name = br.atom_id_2_4c();
                    bool found = false;
                    for (gemmi::Atom const &a : res_ref.atoms) {
                        if (atom_name_to_4c(a.name) == other_atom_name) {
                            found = true;
                            break;
                        }
                    }
                    if (found) {
                        std::string bond_type = br.type();
                        if (bond_type == "single")
                            n_bonds += 1;
                        else if (bond_type == "double")
                            n_bonds += 2;
                        else if (bond_type == "triple")
                            n_bonds += 3;
                    }
                }

                if (br.atom_id_2_4c() == atom_name) {
                    std::string other_atom_name = br.atom_id_1_4c();
                    bool found = false;
                    for (gemmi::Atom const &a : res_ref.atoms) {
                        if (atom_name_to_4c(a.name) == other_atom_name) {
                            found = true;
                            break;
                        }
                    }
                    if (found) {
                        std::string bond_type = br.type();
                        if (bond_type == "single")
                            n_bonds += 1;
                        else if (bond_type == "double")
                            n_bonds += 2;
                        else if (bond_type == "triple")
                            n_bonds += 3;
                    }
                }
            }

            if (n_bonds > 3) { // i.e. 4
                clipper::Coord_orth pt(atom.pos.x, atom.pos.y, atom.pos.z);
                v.push_back(std::make_pair(atom_name, pt));
            }
        }

        return v;
    }
};

} // namespace pli