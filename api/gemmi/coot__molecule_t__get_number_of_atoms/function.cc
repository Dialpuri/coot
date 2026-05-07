#include "function.hh"

namespace coot {
namespace molecule_t {
    unsigned int get_number_of_atoms_gemmi(const gemmi::Structure& st) {
        unsigned int n = 0;
        for (const auto& model : st.models) {
            for (const auto& chain : model.chains) {
                for (const auto& residue : chain.residues) {
                    n += residue.atoms.size();
                }
            }
        }
        return n;
    }
}
}