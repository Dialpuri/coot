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

// atom_name_pair is already defined in get_atom_index_pairs gemmi header,
// but we define atom_name_pair_list_gemmi here since no port exists.

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

std::pair<short int, float> change_by_gemmi(
    int ichi, double diff,
    const std::vector<std::vector<int>>& contact_indices,
    gemmi::Residue& residue)
{
    std::pair<short int, float> p(0, 0.0);

    if (residue.atoms.empty()) {
        return p;
    }

    std::vector<coot::chi_angles::atom_name_pair> atom_name_pairs =
        atom_name_pair_list_gemmi(residue.name);

    if (atom_name_pairs.size() == 1 && atom_name_pairs[0].atom1 == "empty") {
        p.first = 1;
        return p;
    }

    // out of range check
    if (ichi < 1 || ichi > static_cast<int>(atom_name_pairs.size())) {
        p.first = 2;
        p.second = 0.0f;
        return p;
    }

    coot::atom_spec_t dummy_tree_base_atom;
    dummy_tree_base_atom.int_user_data = -999;

    p = change_by_internal_gemmi(ichi, diff, atom_name_pairs,
                                 contact_indices, residue,
                                 dummy_tree_base_atom);

    return p;
}

} // namespace chi_angles
} // namespace coot
