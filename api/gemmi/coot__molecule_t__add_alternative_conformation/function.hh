#pragma once
#include <string>
#include <vector>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>

namespace coot {
namespace molecule_t {

// Helper to parse CID format like "//A/10" into chain_id and residue number
inline std::pair<std::string, int> parse_cid(const std::string& cid) {
    // Format: "//chain/resnum" e.g. "//A/10"
    if (cid.size() < 4 || cid[0] != '/' || cid[1] != '/') {
        return {"", -1};
    }
    size_t slash2 = cid.find('/', 2);
    if (slash2 == std::string::npos) {
        return {"", -1};
    }
    std::string chain_id = cid.substr(2, slash2 - 2);
    int resnum = std::stoi(cid.substr(slash2 + 1));
    return {chain_id, resnum};
}

// Helper to find residue by CID
inline gemmi::Residue* cid_to_residue(gemmi::Model& model, const std::string& cid) {
    auto [chain_id, resnum] = parse_cid(cid);
    if (chain_id.empty() || resnum < 0) {
        return nullptr;
    }
    for (auto& chain : model.chains) {
        if (chain.name == chain_id) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == resnum) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

// Calculate offset vector for alt conf generation
inline gemmi::Vec3 calculate_offset(gemmi::Residue& residue) {
    // Find N and C atoms
    gemmi::Atom* n_at = nullptr;
    gemmi::Atom* c_at = nullptr;
    for (auto& atom : residue.atoms) {
        if (atom.name == " N  ") n_at = &atom;
        if (atom.name == " C  ") c_at = &atom;
    }
    
    if (n_at && c_at) {
        gemmi::Vec3 c_pos(c_at->pos.x, c_at->pos.y, c_at->pos.z);
        gemmi::Vec3 n_pos(n_at->pos.x, n_at->pos.y, n_at->pos.z);
        
        // cn_unit = (c_pos - n_pos).unit()
        gemmi::Vec3 cn_vec = c_pos - n_pos;
        double len = cn_vec.length();
        if (len > 0) {
            cn_vec = cn_vec / len;
        }
        
        // arb_uv = (1,2,3).unit()
        gemmi::Vec3 arb(1, 2, 3);
        double arb_len = arb.length();
        if (arb_len > 0) {
            arb = arb / arb_len;
        }
        
        // cp = cross(cn_unit, arb_uv)
        gemmi::Vec3 cp = cn_vec.cross(arb);
        
        // offset = 0.2 * cp
        return 0.2 * cp;
    }
    
    // Default offset if we can't calculate
    return gemmi::Vec3(0, 0, 0.2);
}

// Main function - ports add_alternative_conformation to gemmi
inline int add_alternative_conformation_gemmi(gemmi::Model& model, const std::string& cid) {
    gemmi::Residue* residue = cid_to_residue(model, cid);
    if (!residue) {
        return -1;  // Residue not found
    }
    
    // Calculate offset
    gemmi::Vec3 offset = calculate_offset(*residue);
    
    std::vector<gemmi::Atom> new_atoms;
    
    for (auto& atom : residue->atoms) {
        if (atom.altloc == '\0') {  // No alt loc currently
            // Create new atom with altloc 'B', shifted by -offset
            gemmi::Atom new_atom = atom;
            new_atom.pos.x -= offset.x;
            new_atom.pos.y -= offset.y;
            new_atom.pos.z -= offset.z;
            new_atom.altloc = 'B';
            new_atom.occ = 0.5;
            new_atoms.push_back(new_atom);
            
            // Modify original atom with altloc 'A', shifted by +offset
            atom.pos.x += offset.x;
            atom.pos.y += offset.y;
            atom.pos.z += offset.z;
            atom.altloc = 'A';
            atom.occ = 0.5;
        }
    }
    
    // Add new atoms to residue
    for (auto& new_atom : new_atoms) {
        residue->atoms.push_back(new_atom);
    }
    
    return 0;  // Success
}

} // namespace molecule_t
} // namespace coot