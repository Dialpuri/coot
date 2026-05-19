#pragma once

#include <vector>
#include <string>
#include <gemmi/model.hpp>

namespace coot {

// Full definition of protein_geometry
class protein_geometry {
public:
    bool have_at_least_minimal_dictionary_for_residue_type_gemmi(
        const std::string& res_name, int imol) const {
        // For standard amino acids and common residues, return true
        // This matches the behavior of the original function
        static const std::vector<std::string> standard_residues = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
            "DA", "DC", "DG", "DT", "A", "C", "G", "T", "U", "DA5", "DC5", "DG5", "DT5",
            "DA3", "DC3", "DG3", "DT3", "A5", "C5", "G5", "T5", "U5", "A3", "C3", "G3", "T3", "U3"
        };
        
        for (const auto& res : standard_residues) {
            if (res_name == res) {
                return true;
            }
        }
        return false;
    }
};

namespace chem_feat_clust {

// chem_feat_solvated_ligand - simplified version with only what's needed
// This mirrors the MMDB version but uses gemmi types
struct chem_feat_solvated_ligand {
    gemmi::Residue* residue;
    // Other members would exist in full version but aren't needed for check_dictionaries
};

// The ported check_dictionaries function
// Checks that dictionaries exist for all ligand residues
inline bool check_dictionaries_gemmi(
    const std::vector<chem_feat_solvated_ligand>& ligands,
    const protein_geometry& geometry_p,
    int imol_fake = 0) {
    
    bool success = true;
    for (unsigned int ilig = 0; ilig < ligands.size(); ilig++) {
        const gemmi::Residue* res = ligands[ilig].residue;
        if (res) {
            std::string res_name = res->name;
            bool have = geometry_p.have_at_least_minimal_dictionary_for_residue_type_gemmi(res_name, imol_fake);
            if (!have) {
                success = false;
                break;
            }
        }
    }
    return success;
}

} // namespace chem_feat_clust
} // namespace coot