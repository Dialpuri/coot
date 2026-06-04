#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <utility>

// Per-atom torsion data extracted from coot::atom_by_torsion_t
struct atom_torsion_data_t {
   std::string atom_name;
   std::string element;
   std::pair<bool, std::string> prior_atom_1;
   std::pair<bool, std::string> prior_atom_2;
   std::pair<bool, std::string> prior_atom_3;
   double bond_length;
   double angle_deg;
   double torsion_deg;
};

namespace coot {

// gemmi port of link_by_torsion_t::make_residue
// Returns nullptr when atom_torsions is empty (equivalent to !filled()).
// Throws std::runtime_error when pos_gemmi cannot find required atoms
// in the base residue (e.g. sugar atoms C1/C2/O5 missing in a protein).
gemmi::Residue* make_residue_gemmi(
    const std::string& new_residue_type,
    int new_res_no,
    const std::vector<atom_torsion_data_t>& atom_torsions,
    const gemmi::Residue& base_residue,
    float b_factor);

} // namespace coot