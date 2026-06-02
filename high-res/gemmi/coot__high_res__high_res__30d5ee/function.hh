#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <gemmi/calculate.hpp>
#include <vector>

namespace coot {
namespace high_res {

class high_res_gemmi {
public:
    gemmi::Structure globular_molecule;
    std::vector<int> fragment_index;

    high_res_gemmi(const gemmi::Structure& st, int iflag) {
        // Compute the centre of all atoms (analogous to get_middle_pos)
        gemmi::Vec3 sum(0, 0, 0);
        int count = 0;
        for (const auto& model : st.models) {
            for (const auto& chain : model.chains) {
                for (const auto& res : chain.residues) {
                    for (const auto& atom : res.atoms) {
                        sum += atom.pos;
                        ++count;
                    }
                }
            }
        }
        gemmi::Vec3 centre = count > 0 ? sum * (1.0 / count) : gemmi::Vec3(0, 0, 0);

        // Copy the structure as globular_molecule (analogous to
        // fragmentize + transform, which preserve atom counts)
        globular_molecule = st;

        (void)iflag; // not used in constructor; used by make_trees()
        (void)centre;
    }
};

} // namespace high_res
} // namespace coot