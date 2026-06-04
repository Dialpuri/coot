#include "function.hh"

void remove_trans_peptide_restraint_gemmi(
    std::vector<simple_restraint_gemmi>& restraints_vec,
    const std::vector<gemmi::Residue*>& atom_residues,
    gemmi::Residue* first,
    gemmi::Residue* second) {

    // Guard against null pointers — match original behavior (no crash)
    if (first == nullptr || second == nullptr) {
        return;
    }

    unsigned int n_rest = static_cast<unsigned int>(restraints_vec.size());
    for (unsigned int i = 0; i < n_rest; i++) {
        simple_restraint_gemmi& restraint = restraints_vec[i];
        if (restraint.restraint_type == TRANS_PEPTIDE_RESTRAINT_VAL) {
            gemmi::Residue* r_11 = atom_residues[restraint.atom_index_1];
            gemmi::Residue* r_12 = atom_residues[restraint.atom_index_2];
            gemmi::Residue* r_21 = atom_residues[restraint.atom_index_3];
            gemmi::Residue* r_22 = atom_residues[restraint.atom_index_4];
            if (r_11 == first) {
                if (r_12 == first) {
                    if (r_21 == second) {
                        if (r_22 == second) {
                            restraint.close();
                        }
                    }
                }
            }
        }
    }
}