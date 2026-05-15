#pragma once

#include <set>
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

/// Inline translation of coot::atom_spec_t using gemmi types.
struct atom_spec_gemmi {
    std::string chain_id;
    int res_no;
    std::string ins_code;
    std::string atom_name;
    std::string alt_conf;

    atom_spec_gemmi() : res_no(0) {}

    atom_spec_gemmi(const gemmi::CRA& cra) : res_no(0) {
        if (cra.chain)
            chain_id = cra.chain->name;
        if (cra.residue) {
            res_no = cra.residue->seqid.num.value;
            ins_code = std::string(1, cra.residue->seqid.icode);
        }
        if (cra.atom) {
            atom_name = cra.atom->name;
            alt_conf = std::string(1, cra.atom->altloc);
        }
    }
};

/// Ported from coot::atom_indices_in_other_molecule.
/// Returns the set of indices into mol_selection where atoms in
/// moving_selection also exist at the same index position, with
/// matching atom name, residue number, and chain ID.
inline std::set<int> atom_indices_in_other_molecule_gemmi(
    const std::vector<gemmi::CRA>& mol_selection,
    const std::vector<gemmi::CRA>& moving_selection) {

    std::set<int> s;

    for (size_t iat = 0; iat < moving_selection.size(); iat++) {
        int idx = static_cast<int>(iat);

        if (idx < static_cast<int>(mol_selection.size())) {
            const auto& mol_cra      = mol_selection[idx];
            const auto& moving_cra   = moving_selection[iat];

            atom_spec_gemmi moving_spec(moving_cra);
            atom_spec_gemmi mol_spec(mol_cra);

            if (moving_spec.atom_name == mol_spec.atom_name &&
                moving_spec.res_no      == mol_spec.res_no &&
                moving_spec.chain_id    == mol_spec.chain_id) {
                s.insert(idx);
            }
        }
    }

    return s;
}

}  // namespace coot