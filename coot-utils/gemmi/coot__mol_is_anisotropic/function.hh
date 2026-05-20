#pragma once
#include <gemmi/model.hpp>
#include <string>

namespace coot {
    // Check if any atom in the structure has anisotropic displacement parameters
    inline bool mol_is_anisotropic_gemmi(const gemmi::Structure& st) {
        bool is_aniso = false;
        if (st.models.empty()) return false;
        
        for (const gemmi::Model& model : st.models) {
            for (const gemmi::Chain& chain : model.chains) {
                for (const gemmi::Residue& residue : chain.residues) {
                    for (const gemmi::Atom& atom : residue.atoms) {
                        if (!atom.aniso.all_zero()) {
                            is_aniso = true;
                            return is_aniso;
                        }
                    }
                }
            }
        }
        return is_aniso;
    }
}