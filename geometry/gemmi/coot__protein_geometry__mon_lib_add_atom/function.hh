#pragma once

#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include "coot/geometry/protein-geometry.hh"

namespace coot {

// Port of coot::protein_geometry::mon_lib_add_atom with MMDB types removed.
// mmdb::realtype → double (same underlying type).
// dict_res_restraints and read_number are private members of protein_geometry,
// so they are passed as explicit parameters.
void mon_lib_add_atom_gemmi(
    std::vector<std::pair<int, dictionary_residue_restraints_t>>& dict_res_restraints,
    int read_number,
    const std::string &comp_id,
    int imol_enc,
    const std::string &atom_id,
    const std::string &atom_id_4c,
    const std::string &type_symbol,
    const std::string &type_energy,
    const std::pair<bool, double> &partial_charge,
    const std::pair<bool, int> &formal_charge,
    dict_atom::aromaticity_t arom_in,
    const std::pair<bool, clipper::Coord_orth> &model_pos,
    const std::pair<bool, clipper::Coord_orth> &model_pos_ideal);

} // namespace coot
