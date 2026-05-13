#include "function.hh"
#include <gemmi/pdb.hpp>
#include <gemmi/chemcomp.hpp>
#include <gemmi/polyheur.hpp>
#include <gemmi/elem.hpp>

namespace coot {

int molecule_t::delete_literal_using_cid_gemmi(const std::string &atom_selection_cids) {
    int status = 0;
    std::vector<gemmi::CRA> atoms_to_be_deleted;
    
    // Split the CID string by "||" delimiter into a vector
    std::vector<std::string> cids;
    std::string remaining = atom_selection_cids;
    size_t pos = 0;
    std::string delimiter = "||";
    while ((pos = remaining.find(delimiter)) != std::string::npos) {
        cids.push_back(remaining.substr(0, pos));
        remaining.erase(0, pos + delimiter.length());
    }
    cids.push_back(remaining);
    
    // Parse each CID and find matching atoms
    for (const auto &cid : cids) {
        // CID format: //chain/resnum/atomname or //chain/resnum-resnum/atomname
        // Remove leading "//" if present
        std::string clean_cid = cid;
        if (clean_cid.size() >= 2 && clean_cid[0] == '/' && clean_cid[1] == '/') {
            clean_cid = clean_cid.substr(2);
        }
        
        // Parse chain/residue/atom format
        std::vector<std::string> parts;
        size_t start = 0;
        for (size_t i = 0; i < clean_cid.size(); ++i) {
            if (clean_cid[i] == '/') {
                parts.push_back(clean_cid.substr(start, i - start));
                start = i + 1;
            }
        }
        if (start < clean_cid.size()) {
            parts.push_back(clean_cid.substr(start));
        }
        
        if (parts.size() < 2) continue;
        
        std::string chain_name = parts[0];
        std::string res_part = parts[1];
        std::string atom_name = parts.size() > 2 ? parts[2] : "";
        
        // Find the chain
        for (auto& model : structure.models) {
            for (auto& chain : model.chains) {
                if (chain.name != chain_name) continue;
                
                // Parse residue range (e.g., "10" or "11-12")
                size_t dash_pos = res_part.find('-');
                int res_start = 0, res_end = 0;
                bool has_range = dash_pos != std::string::npos;
                
                if (has_range) {
                    res_start = std::stoi(res_part.substr(0, dash_pos));
                    res_end = std::stoi(res_part.substr(dash_pos + 1));
                } else {
                    res_start = res_end = std::stoi(res_part);
                }
                
                // Find matching residues and atoms
                for (auto& residue : chain.residues) {
                    if (residue.seqid.num.value >= res_start && 
                        residue.seqid.num.value <= res_end) {
                        
                        if (atom_name.empty()) {
                            // Delete all atoms in residue
                            for (auto& atom : residue.atoms) {
                                atoms_to_be_deleted.push_back({&chain, &residue, &atom});
                            }
                        } else {
                            // Delete specific atom by name
                            for (auto& atom : residue.atoms) {
                                if (atom.name == atom_name) {
                                    atoms_to_be_deleted.push_back({&chain, &residue, &atom});
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    if (!atoms_to_be_deleted.empty()) {
        // Delete atoms by iterating backwards to avoid invalidating iterators
        for (auto it = atoms_to_be_deleted.rbegin(); it != atoms_to_be_deleted.rend(); ++it) {
            const auto& cra = *it;
            if (cra.residue && cra.atom) {
                auto& atoms = cra.residue->atoms;
                for (auto atom_it = atoms.begin(); atom_it != atoms.end(); ++atom_it) {
                    if (&*atom_it == cra.atom) {
                        atoms.erase(atom_it);
                        break;
                    }
                }
            }
        }
        
        status = 1;
        
        // Cleanup - setup entities
        gemmi::setup_entities(structure);
    }
    
    return status;
}

} // namespace coot