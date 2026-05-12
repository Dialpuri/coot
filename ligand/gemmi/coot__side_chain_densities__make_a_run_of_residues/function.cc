#include "function.hh"

namespace coot {

std::vector<gemmi::CRA> side_chain_densities::make_a_run_of_residues_gemmi(const gemmi::Structure& st,
                                                                           const std::string& chain_id,
                                                                           int resno_start, int resno_end) const {
    std::vector<gemmi::CRA> result;
    
    if (st.models.empty()) {
        return result;
    }
    
    const gemmi::Model& model = st.models[0];
    
    for (const gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            for (const gemmi::Residue& residue : chain.residues) {
                int res_no = residue.seqid.num.value;
                if (res_no >= resno_start && res_no <= resno_end) {
                    // Return CRA with chain and residue pointers, atom is nullptr
                    result.push_back({const_cast<gemmi::Chain*>(&chain),
                                     const_cast<gemmi::Residue*>(&residue),
                                     nullptr});
                }
            }
        }
    }
    
    return result;
}

}