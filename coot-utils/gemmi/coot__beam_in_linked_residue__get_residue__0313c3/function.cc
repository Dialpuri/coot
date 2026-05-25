#include "function.hh"

namespace coot {
namespace beam_in_linked_residue {

gemmi::Residue* get_residue_gemmi(const std::string &comp_id, const gemmi::Structure &st) {
    if (st.models.empty())
        return nullptr;
    
    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& residue : chain.residues) {
                if (residue.name == comp_id) {
                    return const_cast<gemmi::Residue*>(&residue);
                }
            }
        }
    }
    
    return nullptr;
}

} // namespace beam_in_linked_residue
} // namespace coot