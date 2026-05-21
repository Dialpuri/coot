#pragma once
#include <utility>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {

struct restraints_container_t {
    enum peptide_order_info_t {
        IS_PEPTIDE,
        IS_NOT_PEPTIDE,
        UNKNOWN
    };
};

inline std::pair<restraints_container_t::peptide_order_info_t, bool>
peptide_C_and_N_are_in_order_p_gemmi(const gemmi::CRA& cra1, const gemmi::CRA& cra2) {
    // Get residues from CRA
    const gemmi::Residue* r1 = cra1.residue;
    const gemmi::Residue* r2 = cra2.residue;
    
    // If residues are null, return UNKNOWN
    if (!r1 || !r2) {
        return std::make_pair(restraints_container_t::UNKNOWN, false);
    }
    
    // Get chain pointers from CRA
    const gemmi::Chain* chain1 = cra1.chain;
    const gemmi::Chain* chain2 = cra2.chain;
    
    // If chains are null, return UNKNOWN
    if (!chain1 || !chain2) {
        return std::make_pair(restraints_container_t::UNKNOWN, false);
    }
    
    // Check if they are on the same chain
    if (chain1 == chain2) {
        // In gemmi, there's no direct index like mmdb::Residue::index
        // We need to find the position in the chain's residue vector
        int idx1 = -1, idx2 = -1;
        for (size_t i = 0; i < chain1->residues.size(); ++i) {
            if (&chain1->residues[i] == r1) idx1 = static_cast<int>(i);
            if (&chain1->residues[i] == r2) idx2 = static_cast<int>(i);
        }
        
        if (idx1 == -1 || idx2 == -1) {
            return std::make_pair(restraints_container_t::UNKNOWN, false);
        }
        
        int serial_delta = idx2 - idx1;
        
        if ((serial_delta == -1) || (serial_delta == 1)) {
            // ok to proceed
        } else {
            return std::make_pair(restraints_container_t::IS_NOT_PEPTIDE, false);
        }
        
        if (serial_delta == 1) {
            // r2 comes after r1 in the chain
            char ins_code_1 = r1->seqid.icode;
            char ins_code_2 = r2->seqid.icode;
            int res_no_delta = r2->seqid.num.value - r1->seqid.num.value;
            
            // Normalize insertion codes: MMDB "" -> gemmi ' ', both mean "no insertion code"
            auto norm_ic = [](char c) { return c == ' ' ? '\0' : c; };
            char norm1 = norm_ic(ins_code_1);
            char norm2 = norm_ic(ins_code_2);
            
            if (norm1 == '\0') {
                if (norm2 == '\0') {
                    if (res_no_delta == 1 || res_no_delta == -1) {
                        return std::make_pair(restraints_container_t::IS_PEPTIDE, false);
                    }
                }
            }
            return std::make_pair(restraints_container_t::UNKNOWN, false);
        } else {
            // serial_delta == -1, r2 comes before r1 in the chain
            char ins_code_1 = r1->seqid.icode;
            char ins_code_2 = r2->seqid.icode;
            int res_no_delta = r2->seqid.num.value - r1->seqid.num.value;
            
            // Normalize insertion codes
            auto norm_ic = [](char c) { return c == ' ' ? '\0' : c; };
            char norm1 = norm_ic(ins_code_1);
            char norm2 = norm_ic(ins_code_2);
            
            if (norm1 == '\0') {
                if (norm2 == '\0') {
                    if (res_no_delta == 1 || res_no_delta == -1) {
                        return std::make_pair(restraints_container_t::IS_PEPTIDE, true);
                    }
                }
            }
            return std::make_pair(restraints_container_t::UNKNOWN, true);
        }
    } else {
        // Different chains - we can't make a decision
        return std::make_pair(restraints_container_t::UNKNOWN, false);
    }
}

} // namespace coot