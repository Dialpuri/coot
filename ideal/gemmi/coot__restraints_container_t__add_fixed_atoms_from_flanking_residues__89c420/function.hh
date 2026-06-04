#pragma once
#include <vector>
#include <set>
#include <cstddef>
#include <gemmi/model.hpp>

namespace coot {

class restraints_container_t_gemmi {
public:
    size_t n_atoms;
    std::vector<gemmi::Residue*> atom_residues;
    std::set<int> fixed_atom_indices;

    void add_fixed_atoms_from_flanking_residues_gemmi(
        bool have_flanking_residue_at_start,
        bool have_flanking_residue_at_end,
        int iselection_start_res,
        int iselection_end_res)
    {
        if (have_flanking_residue_at_start || have_flanking_residue_at_end) {
            for (int iat = 0; iat < static_cast<int>(n_atoms); iat++) {
                gemmi::Residue* res = atom_residues[iat];
                if (have_flanking_residue_at_start) {
                    if (res->seqid.num.value == iselection_start_res) {
                        fixed_atom_indices.insert(iat);
                    }
                }
                if (have_flanking_residue_at_end) {
                    if (res->seqid.num.value == iselection_end_res) {
                        fixed_atom_indices.insert(iat);
                    }
                }
            }
        }
    }
};

} // namespace coot