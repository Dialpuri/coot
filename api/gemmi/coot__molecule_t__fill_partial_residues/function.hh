#pragma once
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <vector>
#include <string>
#include <map>

namespace coot {

struct missing_atom_info_gemmi {
    std::vector<std::string> residues_with_no_dictionary;
    std::vector<gemmi::CRA> residues_with_missing_atoms;
    std::map<gemmi::CRA, std::vector<std::string>> residue_missing_atom_names_map;
    std::vector<std::pair<gemmi::CRA, std::vector<gemmi::Atom*>>> atoms_in_coords_but_not_in_dict;
};

inline missing_atom_info_gemmi missing_atoms(const gemmi::Structure& st,
                                            bool do_missing_hydrogen_atoms_flag,
                                            void* geom_p = nullptr) {
    missing_atom_info_gemmi result;
    (void)do_missing_hydrogen_atoms_flag;
    (void)geom_p;
    return result;
}

inline int fill_partial_residues_gemmi(const gemmi::Structure& st,
                                      void* xmap = nullptr,
                                      void* geom_p = nullptr) {
    auto mai = coot::missing_atoms(st, false, geom_p);
    
    if (!mai.residues_with_missing_atoms.empty()) {
        for (const auto& cra : mai.residues_with_missing_atoms) {
            (void)cra;
        }
    }
    
    return 0;
}

} // namespace coot