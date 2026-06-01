#pragma once

#include <gemmi/model.hpp>

inline constexpr int NO_BOND = 0;

class Bond_lines_container {
public:
    void set_udd_unbonded_gemmi(gemmi::Model const& model) {
        for (gemmi::Chain const& chain : model.chains) {
            for (gemmi::Residue const& res : chain.residues) {
                for (gemmi::Atom const& at : res.atoms) {
                    (void)at;
                }
            }
        }
    }
};