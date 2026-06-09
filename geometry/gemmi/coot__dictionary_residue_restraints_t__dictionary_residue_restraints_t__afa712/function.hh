#pragma once

namespace coot {

struct dictionary_residue_restraints_t {
    bool filled_with_bond_order_data_only_flag;
    bool has_partial_charges_flag;
    int read_number;
    bool nuclear_distances_flag;

    static dictionary_residue_restraints_t dictionary_residue_restraints_t_gemmi(bool constructor_for_srs_restraints) {
        dictionary_residue_restraints_t rest;
        rest.filled_with_bond_order_data_only_flag = true;
        rest.has_partial_charges_flag = false;
        rest.read_number = -1;
        rest.nuclear_distances_flag = false;
        if (constructor_for_srs_restraints) {
            // empty body in original — no-op
        }
        return rest;
    }
};

} // namespace coot
