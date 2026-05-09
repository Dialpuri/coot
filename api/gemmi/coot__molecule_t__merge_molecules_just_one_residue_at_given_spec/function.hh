#pragma once
#include <string>
#include <vector>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {

struct residue_spec_t {
    std::string chain_id;
    int res_no;
    std::string ins_code;
    
    residue_spec_t(const std::string& chain, int res, const std::string& ins) 
        : chain_id(chain), res_no(res), ins_code(ins) {}
    
    residue_spec_t() : chain_id(""), res_no(0), ins_code("") {}
    
    bool empty() const { return chain_id.empty(); }
};

struct molecule_t {
    gemmi::Structure structure;
    gemmi::Model* model() { return &structure.models[0]; }
    
    // Helper to find residue by spec, returns nullptr if not found
    gemmi::Residue* get_residue(const residue_spec_t& spec) {
        for (auto& chain : model()->chains) {
            for (auto& res : chain.residues) {
                if (chain.name == spec.chain_id &&
                    res.seqid.num.value == spec.res_no &&
                    res.seqid.icode == spec.ins_code[0]) {
                    return &res;
                }
            }
        }
        return nullptr;
    }
    
    void make_backup(const std::string& /*unused*/) {}
    
    gemmi::Residue* copy_and_add_residue_to_chain(gemmi::Chain* chain, const gemmi::Residue* src_res) {
        gemmi::Residue new_res = *src_res;
        new_res.atoms.clear();  // Clear atoms, will copy below
        for (const auto& atom : src_res->atoms) {
            gemmi::Atom new_atom = atom;
            new_atom.pos = atom.pos;  // Copy position
            new_res.atoms.push_back(new_atom);
        }
        chain->residues.push_back(new_res);
        return &chain->residues.back();
    }
    
    bool merge_molecules_just_one_residue_at_given_spec_gemmi(
        const gemmi::Structure& molecule_to_add,
        const residue_spec_t& target_spec) {
        
        bool status = false;
        
        if (!target_spec.empty()) {
            gemmi::Residue* residue_p = get_residue(target_spec);
            if (!residue_p) {
                if (model()->chains.empty()) return false;
                
                int n_res = 0;
                for (const auto& model : molecule_to_add.models) {
                    for (const auto& chain : model.chains) {
                        n_res += chain.residues.size();
                    }
                }
                
                if (n_res == 1) {
                    gemmi::Chain* this_chain_p = nullptr;
                    for (auto& chain : model()->chains) {
                        if (chain.name == target_spec.chain_id) {
                            this_chain_p = &chain;
                            break;
                        }
                    }
                    
                    if (!this_chain_p) {
                        gemmi::Chain new_chain;
                        new_chain.name = target_spec.chain_id;
                        model()->chains.push_back(new_chain);
                        this_chain_p = &model()->chains.back();
                    }
                    
                    // Get first residue from molecule_to_add
                    gemmi::Residue* r = nullptr;
                    for (const auto& model : molecule_to_add.models) {
                        for (const auto& chain : model.chains) {
                            if (!chain.residues.empty()) {
                                r = &const_cast<gemmi::Residue&>(chain.residues.front());
                                break;
                            }
                        }
                        if (r) break;
                    }
                    
                    if (r) {
                        make_backup("merge_molecules_just_one_residue_at_given_spec");
                        gemmi::Residue* new_residue_p = copy_and_add_residue_to_chain(this_chain_p, r);
                        new_residue_p->seqid.num.value = target_spec.res_no;
                        if (!target_spec.ins_code.empty()) {
                            new_residue_p->seqid.icode = target_spec.ins_code[0];
                        }
                        status = true;
                    }
                }
            }
        }
        
        if (status) {
            // structure has been modified
        }
        
        return status;
    }
};

} // namespace coot