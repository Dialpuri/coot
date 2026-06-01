#include "function.hh"

namespace coot {

void delete_aniso_records_from_atoms_gemmi(gemmi::Structure &st) {
    for (gemmi::Model &model : st.models) {
        for (gemmi::Chain &chain : model.chains) {
            for (gemmi::Residue &residue : chain.residues) {
                for (gemmi::Atom &atom : residue.atoms) {
                    // Clear anisotropic B-factor data by zeroing all elements
                    atom.aniso = gemmi::SMat33<float>();
                }
            }
        }
    }
}

} // namespace coot