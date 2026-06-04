#include "function.hh"
#include <algorithm>
#include <string>

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

// Trim whitespace helper for gemmi atom name comparison
static std::string trim_atom_name(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// get_atom_index_gemmi — gemmi port of coot::restraints_container_t::get_atom_index
// Searches for an atom by name in the given residue.
// Returns the 0-based index of the atom, or -2 if not found.
// If the residue is not in residues_, falls back to searching in the first residue in residues_.
int restraints_container_t_gemmi::get_atom_index_gemmi(
    const std::string& atom_name_in, gemmi::Residue* residue_p) const {
    
    // Normalize atom name — trim whitespace for comparison since gemmi
    // stores atom names without MMDB-style 4-char padding
    std::string target = trim_atom_name(atom_name_in);
    
    if (!has_valid_residues()) {
        return -2;
    }
    
    // Check if this residue is in the container
    for (size_t i = 0; i < residues_.size(); i++) {
        if (residues_[i].second == residue_p) {
            // Residue found — search its atoms
            for (size_t j = 0; j < residue_p->atoms.size(); j++) {
                if (trim_atom_name(residue_p->atoms[j].name) == target) {
                    return static_cast<int>(j);
                }
            }
            return -2;
        }
    }
    
    // Residue NOT in container — fall back to first residue (same behavior as MMDB original)
    gemmi::Residue* first_res = residues_[0].second;
    for (size_t j = 0; j < first_res->atoms.size(); j++) {
        if (trim_atom_name(first_res->atoms[j].name) == target) {
            return static_cast<int>(j);
        }
    }
    
    return -2;
}

// get_CA_index_gemmi — gemmi port of coot::restraints_container_t::get_CA_index
// Returns the index of the CA atom within the residue.
int restraints_container_t_gemmi::get_CA_index_gemmi(gemmi::Residue* residue_p) const {
    return get_atom_index_gemmi(" CA ", residue_p);
}

} // namespace coot