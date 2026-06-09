#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <utility>
#include <vector>

#include "geometry/residue-and-atom-specs.hh"
#include "geometry/protein-geometry.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__chi_angles__change_by_internal/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_contact_indices_from_restraints__a15507/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_contact_indices_for_PRO_residue/gemmi/function.hh"

namespace coot {

namespace chi_angles {

// Translate atom_name_pair_list — no existing gemmi port
inline std::vector<coot::chi_angles::atom_name_pair>
atom_name_pair_list_gemmi(const std::string& resname) {
    std::vector<coot::chi_angles::atom_name_pair> pairs;
    
    if (resname == "ALA") {
        return pairs; // no chi angles
    } else if (resname == "ARG") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG");
        pairs.emplace_back("CG", "CD");
        pairs.emplace_back("CD", "NE");
        pairs.emplace_back("NE", "CZ");
    } else if (resname == "ASN") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG");
        pairs.emplace_back("CG", "OD1");
    } else if (resname == "ASP") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG");
        pairs.emplace_back("CG", "OD1");
    } else if (resname == "CYS") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "SG");
    } else if (resname == "GLN") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG");
        pairs.emplace_back("CG", "CD");
        pairs.emplace_back("CD", "OE1");
    } else if (resname == "GLU") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG");
        pairs.emplace_back("CG", "CD");
        pairs.emplace_back("CD", "OE1");
    } else if (resname == "GLY") {
        return pairs; // no chi angles
    } else if (resname == "HIS") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG");
        pairs.emplace_back("CG", "ND1");
    } else if (resname == "ILE") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG1");
        pairs.emplace_back("CG1", "CD");
    } else if (resname == "LEU") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG");
        pairs.emplace_back("CG", "CD");
    } else if (resname == "LYS") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG");
        pairs.emplace_back("CG", "CD");
        pairs.emplace_back("CD", "CE");
    } else if (resname == "MET") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG");
        pairs.emplace_back("CG", "SD");
        pairs.emplace_back("SD", "CE");
    } else if (resname == "PHE") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG");
        pairs.emplace_back("CG", "CD");
        pairs.emplace_back("CD", "CE");
    } else if (resname == "PRO") {
        pairs.emplace_back("N", "CA");
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG");
        pairs.emplace_back("CG", "CD");
    } else if (resname == "SER") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "OG");
    } else if (resname == "THR") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "OG1");
    } else if (resname == "TRP") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG");
        pairs.emplace_back("CG", "CD");
        pairs.emplace_back("CD", "NE1");
    } else if (resname == "TYR") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG");
        pairs.emplace_back("CG", "CD");
        pairs.emplace_back("CD", "CE");
    } else if (resname == "VAL") {
        pairs.emplace_back("CA", "CB");
        pairs.emplace_back("CB", "CG1");
        pairs.emplace_back("CB", "CG2");
    } else {
        // Unknown residue
        pairs.emplace_back("empty", "empty");
    }
    
    return pairs;
}

inline std::pair<short int, float>
change_by_gemmi(int ichi, double diff,
                gemmi::Residue* residue,
                coot::protein_geometry* geom_p) {
    
    bool add_reverse_contacts = false;
    
    std::vector<std::vector<int>> contact_indices =
        coot::util::get_contact_indices_from_restraints_gemmi(residue, geom_p, true, add_reverse_contacts);
    
    std::string resname(residue->name);
    
    // PRO variables, used on the outside.
    std::vector<const gemmi::Atom*> ordered_residue_atoms;
    int CA_index = -1; // unfound initially (tree base)
    
    if (resname == "PRO") {
        const gemmi::Atom* n = nullptr;
        const gemmi::Atom* ca = nullptr;
        const gemmi::Atom* cb = nullptr;
        const gemmi::Atom* cg = nullptr;
        const gemmi::Atom* cd = nullptr;
        
        for (const auto& atom : residue->atoms) {
            if (atom.name == " N  ") n = &atom;
            if (atom.name == " CA ") ca = &atom;
            if (atom.name == " CB ") cb = &atom;
            if (atom.name == " CG ") cg = &atom;
            if (atom.name == " CD ") cd = &atom;
        }
        
        if (residue->atoms.size() > 2 && n && ca && cb && cg && cd) {
            // Order: N, CA, CB, CG, CD first, then the rest
            ordered_residue_atoms.reserve(residue->atoms.size());
            ordered_residue_atoms.push_back(n);
            ordered_residue_atoms.push_back(ca);
            ordered_residue_atoms.push_back(cb);
            ordered_residue_atoms.push_back(cg);
            ordered_residue_atoms.push_back(cd);
            
            for (const auto& atom : residue->atoms) {
                if (atom.name != " N  " && atom.name != " CA " &&
                    atom.name != " CB " && atom.name != " CG " &&
                    atom.name != " CD ") {
                    ordered_residue_atoms.push_back(&atom);
                }
            }
        }
        
        if (!ordered_residue_atoms.empty()) {
            contact_indices =
                coot::util::get_contact_indices_for_PRO_residue_gemmi(ordered_residue_atoms, geom_p);
        }
    } // end of specific PRO-logic
    
    std::pair<short int, float> p(0, 0.0);
    
    if (residue->atoms.empty()) {
        p = std::make_pair(2, 0.0);
    } else {
        // Rotatable bonds -> coord indices
        std::vector<coot::chi_angles::atom_name_pair> atom_name_pairs =
            atom_name_pair_list_gemmi(resname);
        
        if (atom_name_pairs.size() == 1) {
            if (atom_name_pairs[0].atom1 == "empty") {
                p.first = 1;
                return p;
            }
        }
        
        // Check ichi bounds
        if (ichi < 1 || ichi > static_cast<int>(atom_name_pairs.size())) {
            p = std::make_pair(2, 0.0);
            return p;
        }
        
        coot::atom_spec_t dummy_tree_base_atom;
        dummy_tree_base_atom.int_user_data = -999; // special "don't use" value
        
        CA_index = 1;
        if (CA_index != -1) {
            dummy_tree_base_atom.int_user_data = CA_index; // ignored!
            dummy_tree_base_atom.atom_name = " N  ";
        }
        
        p = coot::chi_angles::change_by_internal_gemmi(
            ichi, diff, atom_name_pairs,
            contact_indices, *residue,
            dummy_tree_base_atom);
    }
    
    return p;
}

} // namespace chi_angles
} // namespace coot
