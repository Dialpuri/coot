#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <utility>
#include <iostream>

namespace coot {

enum restraint_type_t { BOND_RESTRAINT = 0 };

class restraints_container_t {
public:
    // Pair of (is_fixed, CRA) for each residue
    std::vector<std::pair<bool, gemmi::CRA>> residues_vec;

    struct simple_restraint {
        int type;
        int index1;
        int index2;
        std::vector<bool> fixed_flags;
        double ideal_value;
        double esd;
        double weight;
    };

    std::vector<simple_restraint> restraints_vec;

    restraints_container_t() = default;

    // The ported function
    void make_helix_pseudo_bond_restraints_from_res_vec_gemmi();

    // Helper: returns a (currently empty) fixed_flags vector
    std::vector<bool> make_fixed_flags(int, int) const {
        return std::vector<bool>{};
    }

    // Helper: appends a restraint to restraints_vec
    void add(int type, int index1, int index2, std::vector<bool> fixed_flags,
             double ideal, double esd, double weight) {
        simple_restraint r;
        r.type           = type;
        r.index1         = index1;
        r.index2         = index2;
        r.fixed_flags    = std::move(fixed_flags);
        r.ideal_value    = ideal;
        r.esd            = esd;
        r.weight         = weight;
        restraints_vec.push_back(std::move(r));
    }
};

} // namespace coot

// ---- Inline definition ----

inline void
coot::restraints_container_t::make_helix_pseudo_bond_restraints_from_res_vec_gemmi()
{
    float pseudo_bond_esd = 0.02f;

    for (std::size_t ir = 0; ir < residues_vec.size(); ir++) {
        for (std::size_t jr = 0; jr < residues_vec.size(); jr++) {
            // Same chain? (pointer comparison of the underlying Chain*)
            if (residues_vec[ir].second.chain == residues_vec[jr].second.chain) {
                // At least one is not fixed
                if (residues_vec[ir].first == false || residues_vec[jr].first == false) {
                    int res_no_delta =
                        residues_vec[jr].second.residue->seqid.num.value
                      - residues_vec[ir].second.residue->seqid.num.value;

                    bool jr_is_downstream = false;
                    if (res_no_delta == 3) jr_is_downstream = true;
                    if (res_no_delta == 4) jr_is_downstream = true;

                    if (jr_is_downstream) {
                        // Iterate atoms of residue ir looking for O
                        for (const gemmi::Atom& at_1 :
                             residues_vec[ir].second.residue->atoms) {
                            if (at_1.name != "O") continue;

                            // Iterate atoms of residue jr looking for N or O
                            for (const gemmi::Atom& at_2 :
                                 residues_vec[jr].second.residue->atoms) {
                                // --- O→N restraint ---
                                if (at_2.name == "N") {
                                    if (at_1.altloc == at_2.altloc) {
                                        int index_1 = static_cast<int>(ir);
                                        int index_2 = static_cast<int>(jr);
                                        std::vector<bool> fixed_flags =
                                            make_fixed_flags(index_1, index_2);
                                        double ideal_dist = 2.919;
                                        if (res_no_delta == 3)
                                            ideal_dist = 3.181;
                                        add(BOND_RESTRAINT, index_1, index_2,
                                            fixed_flags, ideal_dist,
                                            pseudo_bond_esd, 1.2);
                                    }
                                }
                                // --- O→O restraint ---
                                if (at_2.name == "O") {
                                    if (at_1.altloc == at_2.altloc) {
                                        int index_1 = static_cast<int>(ir);
                                        int index_2 = static_cast<int>(jr);
                                        std::vector<bool> fixed_flags =
                                            make_fixed_flags(index_1, index_2);
                                        double ideal_dist = 6.16;
                                        if (res_no_delta == 3)
                                            ideal_dist = 4.92;
                                        double O_O_pseudo_bond_esd = 0.07;
                                        add(BOND_RESTRAINT, index_1, index_2,
                                            fixed_flags, ideal_dist,
                                            O_O_pseudo_bond_esd, 1.2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}