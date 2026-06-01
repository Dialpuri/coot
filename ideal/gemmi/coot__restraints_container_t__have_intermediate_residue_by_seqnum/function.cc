#include "function.hh"

namespace coot {

bool have_intermediate_residue_by_seqnum_gemmi(
    gemmi::CRA first,
    gemmi::CRA second,
    const std::vector<std::pair<bool, gemmi::CRA>>& residues_vec) {

    bool r = false;

    if (first.chain && second.chain) {
        if (first.chain->name == second.chain->name) {
            int res_no_1 = first.residue->seqid.num.value;
            int res_no_2 = second.residue->seqid.num.value;
            int res_no_diff = res_no_2 - res_no_1;

            if (res_no_diff != 1) {
                // try to find a residue that has resno more than res_no_1 and
                // less than res_no_2
                for (unsigned int ii = 0; ii < residues_vec.size(); ii++) {
                    int resno_this = residues_vec[ii].second.residue->seqid.num.value;
                    if (resno_this > res_no_1) {
                        if (resno_this < res_no_2) {
                            if (residues_vec[ii].second.chain &&
                                residues_vec[ii].second.chain->name == first.chain->name) {
                                r = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    return r;
}

}  // namespace coot