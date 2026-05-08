#pragma once
#include <vector>
#include <map>
#include <set>
#include <string>
#include <gemmi/pdb.hpp>

namespace coot {

struct molecule_t {
    std::string name;
    int imol_no;
    std::vector<std::pair<bool, gemmi::Residue*>> neighbouring_residues;

    molecule_t(const std::string& name_in, int mol_no_in) : name(name_in), imol_no(mol_no_in) {}
    
    void add_neighbor_residues_for_refinement_help_gemmi(gemmi::Structure& st);
};

}