#include "function.hh"
#include <gemmi/pdb.hpp>
#include <gemmi/math.hpp>

namespace coot {
namespace molecule_t {

std::vector<std::tuple<std::string, int, std::string>>
parse_selection(const std::string& selection_cid, const gemmi::Structure& st) {
    std::vector<std::tuple<std::string, int, std::string>> result;
    
    // Parse selection string like "//A" or "//A/10"
    // Format: //chain/resseq (e.g. //A, //A/10)
    
    if (selection_cid.size() < 2 || selection_cid[0] != '/' || selection_cid[1] != '/')
        return result;
    
    // Extract chain ID (up to '/' or end)
    size_t chain_start = 2;
    size_t chain_end = selection_cid.find('/', chain_start);
    std::string chain_id;
    if (chain_end == std::string::npos) {
        chain_id = selection_cid.substr(chain_start);
    } else {
        chain_id = selection_cid.substr(chain_start, chain_end - chain_start);
    }
    
    // Extract residue sequence number (after second '/')
    int res_seq = -1;
    if (chain_end != std::string::npos && chain_end + 1 < selection_cid.size()) {
        std::string res_str = selection_cid.substr(chain_end + 1);
        try {
            res_seq = std::stoi(res_str);
        } catch (...) {
            res_seq = -1;
        }
    }
    
    // Collect matching atoms
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            if (chain.name != chain_id)
                continue;
            
            for (const auto& res : chain.residues) {
                if (res_seq >= 0 && res.seqid.num.value != res_seq)
                    continue;
                
                for (const auto& atom : res.atoms) {
                    result.push_back(std::make_tuple(chain.name, res.seqid.num.value, atom.name));
                }
            }
        }
    }
    
    return result;
}

int apply_transformation_to_atom_selection_gemmi(
    const std::string& selection_cid,
    int n_atoms,
    const gemmi::Vec3& rotation_centre,
    const gemmi::Transform& rtop,
    gemmi::Structure& st) {
    
    int n_atoms_moved = 0;
    
    // Get the atoms that match the selection
    auto selected_atoms = parse_selection(selection_cid, st);
    
    if (static_cast<int>(selected_atoms.size()) != n_atoms) {
        // Mismatch - return 0 as per original code
        return 0;
    }
    
    // Apply transformation to matching atoms
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                // Check if this residue matches the selection
                bool matches = false;
                for (const auto& sel : selected_atoms) {
                    if (std::get<0>(sel) == chain.name && 
                        std::get<1>(sel) == res.seqid.num.value) {
                        matches = true;
                        break;
                    }
                }
                
                if (matches) {
                    for (auto& atom : res.atoms) {
                        // Apply transformation:
                        // p1 = pt - rotation_centre
                        // p2 = rtop * p1
                        // p3 = p2 - rotation_centre
                        gemmi::Vec3 pt = atom.pos;
                        gemmi::Vec3 p1 = pt - rotation_centre;
                        gemmi::Vec3 p2 = rtop.apply(p1);
                        gemmi::Vec3 p3 = p2 - rotation_centre;
                        atom.pos = gemmi::Position(p3);
                        n_atoms_moved++;
                    }
                }
            }
        }
    }
    
    return n_atoms_moved;
}

} // namespace molecule_t
} // namespace coot