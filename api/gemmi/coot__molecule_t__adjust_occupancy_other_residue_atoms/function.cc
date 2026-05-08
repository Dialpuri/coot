#include "function.hh"
#include <string>
#include <vector>

namespace coot {
namespace molecule_t {

void adjust_occupancy_other_residue_atoms_gemmi(
    gemmi::Atom* atom,
    gemmi::Residue* residue,
    bool force_sum_1_flag) {

    if (!force_sum_1_flag) {
        std::string new_atom_name = atom->name;
        char new_atom_altloc = atom->altloc;
        std::vector<gemmi::Atom*> same_name_atoms;
        float sum_occ = 0;
        
        for (gemmi::Atom& other_atom : residue->atoms) {
            if (other_atom.name == new_atom_name) {
                if (other_atom.altloc != new_atom_altloc) {
                    same_name_atoms.push_back(&other_atom);
                    sum_occ += other_atom.occ;
                }
            }
        }

        if (sum_occ > 0.01f && same_name_atoms.size() > 0) {
            float other_atom_occ_sum = 0.0f;
            for (size_t i = 0; i < same_name_atoms.size(); i++)
                other_atom_occ_sum += same_name_atoms[i]->occ;

            float remainder = 1.0f - atom->occ;
            float f = remainder / other_atom_occ_sum;
            for (size_t i = 0; i < same_name_atoms.size(); i++) {
                same_name_atoms[i]->occ *= f;
            }
        }
    }
}

}
}