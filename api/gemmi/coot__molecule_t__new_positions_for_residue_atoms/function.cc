#include "function.hh"
#include <gemmi/pdb.hpp>
#include <string>
#include <vector>
#include <map>

namespace coot {
namespace molecule_t {

// Helper: parse CID string like "//A/10" to get chain name and residue number
static std::pair<std::string, int> parse_cid(const std::string& cid) {
    // CID format: //chain/resnum (e.g., "//A/10")
    // The format starts with // so we need to handle that
    // first_slash is at position 0, second_slash is at position 1
    // chain_name starts after second_slash
    // residue number is after the third slash
    
    size_t first_slash = cid.find('/');
    size_t second_slash = cid.find('/', first_slash + 1);
    
    if (second_slash == std::string::npos || first_slash == std::string::npos) {
        return {"", -1};
    }
    
    // Find the third slash (after chain name)
    size_t third_slash = cid.find('/', second_slash + 1);
    
    std::string chain_name;
    if (third_slash != std::string::npos) {
        // Extract chain name between second_slash and third_slash
        chain_name = cid.substr(second_slash + 1, third_slash - second_slash - 1);
    } else {
        // No third slash found, try to parse as much as possible after second_slash
        chain_name = cid.substr(second_slash + 1);
    }
    
    std::string resnum_str;
    if (third_slash != std::string::npos) {
        resnum_str = cid.substr(third_slash + 1);
    } else {
        // No third slash - extract the number at the end
        resnum_str = cid.substr(second_slash + 1);
    }
    
    // Trim whitespace from residue number string
    size_t start_pos = resnum_str.find_first_not_of(" \t\r\n");
    if (start_pos == std::string::npos) {
        return {"", -1};
    }
    size_t end_pos = resnum_str.find_last_not_of(" \t\r\n");
    resnum_str = resnum_str.substr(start_pos, end_pos - start_pos + 1);
    
    int resnum = std::stoi(resnum_str);
    
    return {chain_name, resnum};
}

int new_positions_for_residue_atoms_gemmi(
    gemmi::Model& model,
    const std::string& residue_cid,
    const std::vector<coot::api::moved_atom_t>& moved_atoms) {
    
    try {
        // Parse the CID to get chain name and residue number
        auto [chain_name, resnum] = parse_cid(residue_cid);
        
        if (chain_name.empty() || resnum < 0) {
            return -1;
        }
        
        // Find the chain
        gemmi::Chain* chain = nullptr;
        for (auto& c : model.chains) {
            if (c.name == chain_name) {
                chain = &c;
                break;
            }
        }
        
        if (!chain) {
            return -1;
        }
        
        // Find the residue
        gemmi::Residue* res = nullptr;
        for (auto& r : chain->residues) {
            if (r.seqid.num.value == resnum) {
                res = &r;
                break;
            }
        }
        
        if (!res) {
            return -1;
        }
        
        // Create a map from atom name to index for lookup
        // Need to handle potential whitespace in atom names (gemmi keeps names as-is from PDB)
        std::map<std::string, int> atom_name_to_idx;
        for (size_t i = 0; i < res->atoms.size(); ++i) {
            // Trim leading/trailing spaces from atom name for lookup
            std::string trimmed_name = res->atoms[i].name;
            size_t start = trimmed_name.find_first_not_of(" \t");
            size_t end = trimmed_name.find_last_not_of(" \t");
            if (start != std::string::npos) {
                trimmed_name = trimmed_name.substr(start, end - start + 1);
            }
            atom_name_to_idx[trimmed_name] = i;
        }
        
        // Update positions for each moved atom
        int processed = 0;
        for (const auto& moved : moved_atoms) {
            // Trim the moved atom name as well
            std::string trimmed_moved_name = moved.atom_name;
            size_t start = trimmed_moved_name.find_first_not_of(" \t");
            size_t end = trimmed_moved_name.find_last_not_of(" \t");
            if (start != std::string::npos) {
                trimmed_moved_name = trimmed_moved_name.substr(start, end - start + 1);
            }
            
            auto it = atom_name_to_idx.find(trimmed_moved_name);
            if (it != atom_name_to_idx.end()) {
                // Update the position
                res->atoms[it->second].pos.x = moved.x;
                res->atoms[it->second].pos.y = moved.y;
                res->atoms[it->second].pos.z = moved.z;
                processed++;
            }
        }
        
        return processed;
    } catch (const std::exception& e) {
        return -1;
    }
}

} // namespace molecule_t
} // namespace coot