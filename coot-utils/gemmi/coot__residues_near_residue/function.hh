#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/neighbor.hpp>
#include <set>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

// Parse CID string like "//A/10" to extract chain_id, res_no, ins_code
inline residue_spec_t cid_to_residue_spec(const std::string &cid) {
    // Format: "//A/10" or "//A/10A" etc.
    // Remove leading "//"
    std::string s = cid;
    if (s.size() >= 2 && s[0] == '/' && s[1] == '/') {
        s = s.substr(2);
    }
    
    // Find first '/' for chain ID
    size_t pos1 = s.find('/');
    if (pos1 == std::string::npos) {
        return residue_spec_t(); // invalid
    }
    
    std::string chain_id = s.substr(0, pos1);
    
    // Get remainder after chain ID
    std::string remainder = s.substr(pos1 + 1);
    
    // Parse residue number and optional insertion code
    int res_no = 0;
    std::string ins_code_str = "";
    
    // Find where digits end
    size_t i = 0;
    while (i < remainder.size() && (isdigit(remainder[i]) || remainder[i] == '-' || remainder[i] == '+')) {
        i++;
    }
    
    res_no = std::stoi(remainder.substr(0, i));
    if (i < remainder.size()) {
        ins_code_str = remainder.substr(i);
    }
    
    // Default to space if no insertion code (empty string)
    if (ins_code_str.empty()) {
        ins_code_str = " ";
    }
    
    return residue_spec_t(chain_id, res_no, ins_code_str);
}

// Ported version of residues_near_residue that takes gemmi::Structure instead of mmdb::Manager
inline std::vector<coot::residue_spec_t>
residues_near_residue_gemmi(const coot::residue_spec_t &rs,
                            gemmi::Structure &st,
                            float radius) {
    std::vector<coot::residue_spec_t> result;
    
    // Find the reference residue
    gemmi::Residue *ref_residue = nullptr;
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                if (chain.name == rs.chain_id &&
                    res.seqid.num.value == rs.res_no &&
                    res.seqid.icode == rs.ins_code[0]) {
                    ref_residue = &res;
                    break;
                }
            }
            if (ref_residue) break;
        }
        if (ref_residue) break;
    }
    
    if (!ref_residue) {
        // Residue not found, return empty vector
        return result;
    }
    
    // Get all atoms from the reference residue
    std::vector<gemmi::Atom*> ref_atoms;
    for (gemmi::Atom& atom : ref_residue->atoms) {
        ref_atoms.push_back(&atom);
    }
    
    if (ref_atoms.empty()) {
        return result;
    }
    
    // Build neighbor search structure
    gemmi::NeighborSearch ns(st.models[0], st.cell, radius);
    ns.populate(/*include_h=*/false);
    
    // Find all contacts from reference residue atoms
    std::set<std::string> nearby_residue_keys;  // Use chain_name:res_seq:icode as key
    
    // Manually iterate through reference residue atoms and find contacts
    for (gemmi::Atom* ref_atom : ref_atoms) {
        std::vector<gemmi::NeighborSearch::Mark*> hits =
            ns.find_atoms(gemmi::Position(ref_atom->pos), '\0', 0.0, radius);
        
        for (auto* mark : hits) {
            gemmi::CRA cra = mark->to_cra(st.models[0]);
            if (cra.residue && cra.residue != ref_residue) {
                // Create a unique key for this residue
                std::string key = cra.chain->name + ":" +
                                std::to_string(cra.residue->seqid.num.value) + ":" +
                                std::string(1, cra.residue->seqid.icode);
                nearby_residue_keys.insert(key);
            }
        }
    }
    
    // Convert keys back to residue_spec_t
    for (const auto& key : nearby_residue_keys) {
        // Parse the key back
        size_t pos1 = key.find(':');
        size_t pos2 = key.find(':', pos1 + 1);
        
        if (pos1 != std::string::npos && pos2 != std::string::npos) {
            std::string chain_id = key.substr(0, pos1);
            int res_no = std::stoi(key.substr(pos1 + 1, pos2 - pos1 - 1));
            std::string ins_code = key.substr(pos2 + 1);
            
            if (ins_code.empty()) {
                ins_code = " ";
            }
            
            result.push_back(coot::residue_spec_t(chain_id, res_no, ins_code));
        }
    }
    
    return result;
}

} // namespace coot