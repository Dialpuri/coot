#include "function.hh"

namespace coot {

restraints_container_t_gemmi::restraints_container_t_gemmi(
    const std::vector<std::pair<bool, gemmi::Residue*>>& residues,
    const std::string& geom_type,
    gemmi::Structure* st,
    const clipper::Xmap<float>* map_p_in) {
    residues_.clear();
    for (size_t i = 0; i < residues.size(); i++) {
        if (residues[i].second) {
            residues_.push_back(residues[i]);
        }
    }

    std::sort(residues_.begin(), residues_.end(), 
              [](const std::pair<bool, gemmi::Residue*>& a, 
                 const std::pair<bool, gemmi::Residue*>& b) {
                  return a.second->seqid.num.value < b.second->seqid.num.value;
              });
}

bool restraints_container_t_gemmi::has_valid_residues() const {
    for (size_t i = 0; i < residues_.size(); i++) {
        if (residues_[i].second) {
            return true;
        }
    }
    return false;
}

} // namespace coot