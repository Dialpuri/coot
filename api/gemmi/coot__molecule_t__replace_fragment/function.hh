#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {

struct atom_selection_container_t {
    gemmi::Model* model;
    std::vector<std::tuple<gemmi::Atom*, std::string, int, char>> atoms_with_context;
    // tuple: (atom, chain_name, seq_num, ins_code)
    int UDDOldAtomIndexHandle;
    
    atom_selection_container_t() : model(nullptr), UDDOldAtomIndexHandle(-1) {}
};

struct residue_spec_t {
    std::string chain_name;
    int seq_num;
    char ins_code;
    
    residue_spec_t(const std::string& chain, int seq, char ins) 
        : chain_name(chain), seq_num(seq), ins_code(ins) {}
};

// Helper function to find residue in gemmi model
inline gemmi::Residue* find_residue(gemmi::Model& model, const residue_spec_t& spec) {
    for (auto& chain : model.chains) {
        if (chain.name == spec.chain_name) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == spec.seq_num && res.seqid.icode == spec.ins_code) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

// Helper function to find chain in gemmi model
inline gemmi::Chain* find_chain(gemmi::Model& model, const std::string& chain_id) {
    for (auto& chain : model.chains) {
        if (chain.name == chain_id) {
            return &chain;
        }
    }
    return nullptr;
}

// Ported replace_fragment function - replaces fragment atoms in the model
// Returns 1 on success
int replace_fragment_gemmi(gemmi::Model& model, atom_selection_container_t asc);

// Helper to get chain by ID (from the original lambda)
inline gemmi::Chain* get_chain(gemmi::Model& model, const std::string& chain_id) {
    return find_chain(model, chain_id);
}

// Helper to get residue (from original function)
inline gemmi::Residue* get_residue(gemmi::Model& model, const residue_spec_t& spec) {
    return find_residue(model, spec);
}

// Helper to find serial number for insert (simplified version)
inline int find_serial_number_for_insert_gemmi(int res_no, const std::string& ins_code, const std::string& chain_id) {
    // Simplified - in real coot this uses internal ordering
    return -1;  // Return -1 to mean "not found", will add at end
}

// Helper to make backup (stubbed - no-op in ported version)
inline void make_backup_gemmi(const std::string& /*operation*/) {
    // In real coot this maintains undo history
}

} // namespace coot