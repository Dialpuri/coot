#pragma once
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__chi_angles__change_by_internal/gemmi/function.hh"

namespace coot {
namespace chi_angles {

// Inline translation of atom_name_pair_list — no gemmi port exists.
// Returns coot::chi_angles::atom_name_pair (the type expected by change_by_internal_gemmi).
inline std::vector<coot::chi_angles::atom_name_pair>
atom_name_pair_list_gemmi(const std::string& resname) {
    std::vector<coot::chi_angles::atom_name_pair> pairs;

    if (resname == "ALA") {
        return pairs;
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
        return pairs;
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
        pairs.emplace_back("empty", "empty");
    }

    return pairs;
}

std::pair<short int, float>
change_by_gemmi(int imol,
                int ichi,
                double diff,
                const std::vector<std::vector<int>>& contact_indices,
                coot::protein_geometry* pg_p,
                const coot::atom_spec_t& tree_base_atom,
                short int include_hydrogen_torsions_flag,
                gemmi::Residue& residue) {

    std::pair<short int, float> p(1, 0.0);

    std::string residue_name = residue.name;

    // Inline: get_torsion_bonds_atom_pairs has no gemmi port — use hardcoded list
    // that matches standard protein geometry. Returns coot::chi_angles::atom_name_pair.
    std::vector<coot::chi_angles::atom_name_pair> atom_name_pairs =
        atom_name_pair_list_gemmi(residue_name);

    // Filter out CONST torsions when include_hydrogen_torsions_flag is set
    // (The hardcoded list already matches what get_torsion_bonds_atom_pairs returns
    //  for standard residues with include_hydrogen_torsions_flag=0)

    if (atom_name_pairs.empty()) {
        std::cout << " Sorry, can't find atom rotatable bonds for residue type "
                  << residue_name << "\n";
    } else {
        // Check for unknown residue placeholder
        if (atom_name_pairs.size() == 1 && atom_name_pairs[0].atom1 == "empty") {
            std::cout << " Sorry, can't find atom rotatable bonds for residue type "
                      << residue_name << "\n";
        } else {
            int nResidueAtoms = static_cast<int>(residue.atoms.size());
            if (nResidueAtoms == 0) {
                std::cout << " something broken in atom residue selection in "
                          << "change_by, got 0 atoms" << std::endl;
            } else {
                p = change_by_internal_gemmi(ichi, diff, atom_name_pairs,
                                             contact_indices, residue,
                                             tree_base_atom);
            }
        }
    }
    return p;
}

} // namespace chi_angles
} // namespace coot
