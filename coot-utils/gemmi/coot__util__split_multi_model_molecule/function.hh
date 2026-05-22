#pragma once
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/mmread.hpp>

namespace coot {
namespace util {

// Split a multi-model structure into individual structures, each containing one model
std::vector<gemmi::Structure>
split_multi_model_molecule_gemmi(const gemmi::Structure& st) {
    std::vector<gemmi::Structure> result;
    
    for (const gemmi::Model& model : st.models) {
        gemmi::Structure new_st = st.empty_copy();
        gemmi::Model new_model = model.empty_copy();
        // Copy all atoms from the original model to the new model
        for (const gemmi::Chain& chain : model.chains) {
            gemmi::Chain new_chain = chain.empty_copy();
            for (const gemmi::Residue& residue : chain.residues) {
                gemmi::Residue new_residue = residue.empty_copy();
                for (const gemmi::Atom& atom : residue.atoms) {
                    new_residue.atoms.push_back(atom);
                }
                new_chain.residues.push_back(new_residue);
            }
            new_model.chains.push_back(new_chain);
        }
        new_st.models.push_back(new_model);
        result.push_back(new_st);
    }
    
    return result;
}

} // namespace util
} // namespace coot