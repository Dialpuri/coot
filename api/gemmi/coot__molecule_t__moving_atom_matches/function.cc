#include "function.hh"
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {

bool moving_atom_matches_gemmi(const std::vector<gemmi::CRA>& selection,
                               size_t ref_idx,
                               const gemmi::CRA& moving_cra) {
    bool matches = false;
    
    if (selection.size() > 0) {
        if (ref_idx >= selection.size()) {
            return false;
        } else {
            const gemmi::CRA& ref_cra = selection[ref_idx];
            
            // Extract moving atom properties
            std::string atom_name_mov = moving_cra.atom->name;
            char ins_code_mov = moving_cra.residue->seqid.icode;
            char alt_conf_mov = moving_cra.atom->altloc;
            std::string chain_id_mov = moving_cra.chain ? moving_cra.chain->name : "";
            int resno_mov = moving_cra.residue->seqid.num.value;
            
            // Extract reference atom properties
            std::string atom_name_ref = ref_cra.atom->name;
            char ins_code_ref = ref_cra.residue->seqid.icode;
            char alt_conf_ref = ref_cra.atom->altloc;
            std::string chain_id_ref = ref_cra.chain ? ref_cra.chain->name : "";
            int resno_ref = ref_cra.residue->seqid.num.value;
            
            if (atom_name_ref == atom_name_mov) {
                if (ins_code_ref == ins_code_mov) {
                    if (resno_ref == resno_mov) {
                        if (alt_conf_ref == alt_conf_mov) {
                            if (chain_id_mov == chain_id_ref) {
                                matches = true;
                            }
                        }
                    }
                }
            }
        }
    }
    return matches;
}

} // namespace molecule_t
} // namespace coot