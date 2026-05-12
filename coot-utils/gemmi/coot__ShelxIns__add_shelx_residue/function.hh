#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include "coot/coot-utils/coot-shelx.hh"

namespace coot {

// Port of coot::ShelxIns::add_shelx_residue to gemmi
// Returns a newly allocated gemmi::Residue*
inline gemmi::Residue*
add_shelx_residue_gemmi(const ShelxIns& /*si*/, std::vector<gemmi::Atom>& atom_vector,
                        const std::string& current_res_name, int current_res_no) {
    gemmi::Residue* residue = new gemmi::Residue();
    residue->name = current_res_name;
    residue->seqid.num.value = current_res_no;
    
    // Determine if residue is standard (for Het flag)
    auto is_standard_residue_name = [](const std::string& res_name) {
        // Standard amino acids
        if (res_name == "ALA" || res_name == "ARG" || res_name == "ASN" ||
            res_name == "ASP" || res_name == "CYS" || res_name == "GLN" ||
            res_name == "GLU" || res_name == "GLY" || res_name == "HIS" ||
            res_name == "ILE" || res_name == "LEU" || res_name == "LYS" ||
            res_name == "MET" || res_name == "MSE" || res_name == "PHE" ||
            res_name == "PRO" || res_name == "SER" || res_name == "THR" ||
            res_name == "TRP" || res_name == "TYR" || res_name == "VAL") {
            return true;
        }
        // Standard nucleotides and variants
        if (res_name == "G" || res_name == "A" || res_name == "T" ||
            res_name == "C" || res_name == "U" || res_name == "DG" ||
            res_name == "DA" || res_name == "DT" || res_name == "DC" ||
            res_name == "DU" || res_name == "GR" || res_name == "AR" ||
            res_name == "UR" || res_name == "TR" || res_name == "Gr" ||
            res_name == "Ar" || res_name == "Ur" || res_name == "Tr" ||
            res_name == "Gd" || res_name == "Ad" || res_name == "Ud" ||
            res_name == "Td" || res_name == "G " || res_name == "A " ||
            res_name == "T " || res_name == "C " || res_name == "U " ||
            res_name == "DG " || res_name == "DA " || res_name == "DT " ||
            res_name == "DC ") {
            return true;
        }
        return false;
    };
    
    bool srn = is_standard_residue_name(current_res_name);
    
    for (gemmi::Atom& atom : atom_vector) {
        if (!srn) {
            atom.occ = 0; // Het=1 in original MMDB
        }
        // Atom is already in the vector (no AddAtom needed)
    }
    
    return residue;
}

} // namespace coot