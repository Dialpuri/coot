#pragma once
#include <vector>
#include <utility>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <clipper/core/xmap.h>

namespace coot {

class restraints_container_t_gemmi {
public:
    restraints_container_t_gemmi(const std::vector<std::pair<bool, gemmi::Residue*>>& residues,
                                 const std::string& geom_type,
                                 gemmi::Structure* st,
                                 const clipper::Xmap<float>* map_p_in);
    
    size_t get_residues_vec_size() const { return residues_.size(); }
    bool has_valid_residues() const;

    // get_CA_index_gemmi — gemmi port of coot::restraints_container_t::get_CA_index
    // Returns the index of the CA atom within the residue's atom list,
    // or -2 if not found.
    int get_CA_index_gemmi(gemmi::Residue* residue_p) const;

    // get_atom_index_gemmi — helper, gemmi port of coot::restraints_container_t::get_atom_index
    int get_atom_index_gemmi(const std::string& atom_name_in, gemmi::Residue* residue_p) const;

private:
    std::vector<std::pair<bool, gemmi::Residue*>> residues_;
};

} // namespace coot