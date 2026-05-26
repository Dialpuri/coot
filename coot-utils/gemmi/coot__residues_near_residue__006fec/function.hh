#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/neighbor.hpp>
#include <set>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

inline std::vector<coot::residue_spec_t>
residues_near_residue_gemmi(const coot::residue_spec_t &rs,
                            gemmi::Structure &st,
                            float radius) {
    std::vector<coot::residue_spec_t> result;
    
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
        return result;
    }
    
    std::vector<gemmi::Atom*> ref_atoms;
    for (gemmi::Atom& atom : ref_residue->atoms) {
        ref_atoms.push_back(&atom);
    }
    
    if (ref_atoms.empty()) {
        return result;
    }
    
    gemmi::NeighborSearch ns(st.models[0], st.cell, radius);
    ns.populate(false);
    
    std::set<std::string> nearby_residue_keys;
    
    for (gemmi::Atom* ref_atom : ref_atoms) {
        std::vector<gemmi::NeighborSearch::Mark*> hits =
            ns.find_atoms(gemmi::Position(ref_atom->pos), '\0', 0.0, radius);
        
        for (auto* mark : hits) {
            gemmi::CRA cra = mark->to_cra(st.models[0]);
            if (cra.residue && cra.residue != ref_residue) {
                std::string key = cra.chain->name + ":" +
                                std::to_string(cra.residue->seqid.num.value) + ":" +
                                std::string(1, cra.residue->seqid.icode);
                nearby_residue_keys.insert(key);
            }
        }
    }
    
    for (const auto& key : nearby_residue_keys) {
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