#pragma once
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <string>
#include <algorithm>

namespace coot {

// Port of coot::unshelx - converts a ShelX-format molecule to standard format
// Splits a single chain into multiple chains based on residue numbering gaps
// Returns a new gemmi::Structure with potentially multiple chains
// Returns nullptr if input is invalid (null, wrong number of chains, etc.)
inline gemmi::Structure *unshelx_gemmi(const gemmi::Structure *shelx_st) {
    if (!shelx_st) {
        return nullptr;
    }
    
    // Only process first model
    if (shelx_st->models.empty()) {
        return nullptr;
    }
    
    const gemmi::Model& shelx_model = shelx_st->models[0];
    
    // Check for exactly one chain
    int nchains = 0;
    for (const auto& chain : shelx_model.chains) {
        nchains++;
    }
    
    if (nchains != 1) {
        return nullptr;
    }
    
    gemmi::Structure *mol = new gemmi::Structure();
    mol->name = shelx_st->name;
    mol->connections = shelx_st->connections;
    mol->cell = shelx_st->cell;
    mol->spacegroup_hm = shelx_st->spacegroup_hm;
    
    gemmi::Model model;
    int skip_chain_step = 21;
    std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
    int r_index = 0;
    
    const gemmi::Chain& shelx_chain = shelx_model.chains[0];
    bool need_new_chain = true;
    int ires_prev = -1000;
    
    for (const gemmi::Residue& shelx_res : shelx_chain.residues) {
        int resno = shelx_res.seqid.num.value;
        
        if (resno > (ires_prev + skip_chain_step)) {
            need_new_chain = true;
        }
        
        if (need_new_chain) {
            gemmi::Chain chain;
            std::string new_chain_id = r.substr(r_index, 1);
            r_index++;
            chain.name = new_chain_id;
            model.chains.push_back(chain);
            need_new_chain = false;
        }
        
        // Add a copy of the residue
        gemmi::Residue copy_res = shelx_res;
        model.chains.back().residues.push_back(copy_res);
        
        ires_prev = resno;
    }
    
    mol->models.push_back(model);
    
    // Fix residue indices
    for (gemmi::Chain& chain : mol->models[0].chains) {
        for (size_t ires = 0; ires < chain.residues.size(); ++ires) {
            chain.residues[ires].seqid.num.value = static_cast<int>(ires + 1);
        }
    }
    
    // Remove empty chains
    mol->models[0].chains.erase(
        std::remove_if(mol->models[0].chains.begin(), mol->models[0].chains.end(),
                      [](const gemmi::Chain& c) { return c.residues.empty(); }),
        mol->models[0].chains.end()
    );
    
    return mol;
}

}