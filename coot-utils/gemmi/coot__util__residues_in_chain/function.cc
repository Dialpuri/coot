#include "function.hh"
#include <gemmi/model.hpp>

namespace coot { namespace util {

std::vector<gemmi::Residue*> residues_in_chain_gemmi(gemmi::Model& model, const std::string& chain_id) {
    std::vector<gemmi::Residue*> v;
    
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            for (gemmi::Residue& residue : chain.residues) {
                v.push_back(&residue);
            }
        }
    }
    
    return v;
}

}} // namespace coot::util