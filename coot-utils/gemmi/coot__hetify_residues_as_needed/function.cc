#include "function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__hetify_residue_atoms_as_needed/gemmi/function.hh"

namespace coot {

int hetify_residues_as_needed_gemmi(gemmi::Model* model) {
    int r = 0;
    if (model) {
        for (gemmi::Chain& chain : model->chains) {
            for (gemmi::Residue& res : chain.residues) {
                coot::hetify_residue_atoms_as_needed_gemmi(&res);
            }
        }
    }
    return r;
}

} // namespace coot