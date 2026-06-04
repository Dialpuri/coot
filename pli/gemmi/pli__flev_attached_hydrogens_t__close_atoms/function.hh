#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <vector>

namespace pli {

class flev_attached_hydrogens_t {
public:
    std::vector<gemmi::CRA> close_atoms_gemmi(
        const gemmi::Vec3 &pt,
        const std::vector<gemmi::Residue*> &env_residues,
        gemmi::Model &model) const {

        std::vector<gemmi::CRA> v;
        double dist_crit = 6.0;
        double dist_crit_squared = dist_crit * dist_crit;

        for (gemmi::Chain &chain : model.chains) {
            for (gemmi::Residue &res : chain.residues) {
                bool is_env = false;
                for (gemmi::Residue *r : env_residues) {
                    if (r == &res) {
                        is_env = true;
                        break;
                    }
                }
                if (!is_env) continue;
                if (res.name == "HOH") continue;

                for (gemmi::Atom &atom : res.atoms) {
                    gemmi::Vec3 diff = pt - atom.pos;
                    double d_squared = diff.length_sq();
                    if (d_squared < dist_crit_squared) {
                        v.push_back(gemmi::CRA{&chain, &res, &atom});
                    }
                }
            }
        }
        return v;
    }
};

} // namespace pli