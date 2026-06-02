#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

void update_chem_comp_atoms_from_residue_gemmi(
    const gemmi::Residue& residue,
    dictionary_residue_restraints_t* restraints) {

    for (const gemmi::Atom& at : residue.atoms) {
        std::string atom_name(at.name);

        std::vector<dict_atom>& atom_info = restraints->atom_info;
        for (unsigned int jat = 0; jat < atom_info.size(); jat++) {
            dict_atom& da = atom_info[jat];
            if (da.atom_id_4c == atom_name) {
                clipper::Coord_orth c = coot::co_gemmi(&at);
                da.model_Cartn.second = c;
            }
        }
    }
}

} // namespace coot