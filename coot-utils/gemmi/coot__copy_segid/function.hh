#pragma once
#include <gemmi/model.hpp>
#include <gemmi/polyheur.hpp>
#include <string>
#include <stdexcept>

namespace coot {

std::string residue_atoms_segid_gemmi(const gemmi::Residue& residue) {
    std::vector<std::string> seg_ids;
    
    for (const gemmi::Atom& atom : residue.atoms) {
        std::string seg_id = residue.segment;
        if (seg_ids.empty()) {
            seg_ids.push_back(seg_id);
        } else {
            bool found = false;
            for (const auto& s : seg_ids) {
                if (s == seg_id) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::string mess = "No consistent segids for residue ";
                mess += std::to_string(residue.seqid.num.value);
                throw std::runtime_error(mess);
            }
        }
    }
    
    if (seg_ids.empty()) {
        std::string mess = "No segids for residue ";
        mess += std::to_string(residue.seqid.num.value);
        throw std::runtime_error(mess);
    }
    
    return seg_ids[0];
}

// Copy segid from provider residue to receiver residue
// Returns true on success, false on error
inline bool copy_segid_gemmi(const gemmi::Residue& provider, gemmi::Residue& receiver) {
    try {
        std::string s = residue_atoms_segid_gemmi(provider);
        // In gemmi, segid is stored in ResidueId::segment
        receiver.segment = s.substr(0, 4);
        return true;
    } catch (const std::runtime_error& mess) {
        return false;
    }
}

} // namespace coot