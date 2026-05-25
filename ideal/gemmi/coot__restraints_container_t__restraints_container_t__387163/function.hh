#pragma once
#include <vector>
#include <utility>
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
    
private:
    std::vector<std::pair<bool, gemmi::Residue*>> residues_;
};

} // namespace coot