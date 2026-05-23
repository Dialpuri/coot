#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {
namespace reduce {

bool
hack_ss_bond_test_gemmi(const gemmi::Residue* residue_p, const gemmi::Model& model) {
    bool status = false;

    const gemmi::Atom* CYS_1_SG = nullptr;
    for (const auto& at : residue_p->atoms) {
        if (at.name == "SG") {
            CYS_1_SG = &at;
            break;
        }
    }
    if (!CYS_1_SG) return false;

    clipper::Coord_orth pt_1 = co_gemmi(CYS_1_SG);
    
    for (const auto& chain : model.chains) {
        for (const auto& res : chain.residues) {
            if (res.name == "CYS" && &res != residue_p) {
                for (const auto& at : res.atoms) {
                    if (at.name == "SG") {
                        clipper::Coord_orth pt_2 = co_gemmi(&at);
                        double dd = (pt_2 - pt_1).lengthsq();
                        if (dd < 3.0 * 3.0) {
                            status = true;
                            return status;
                        }
                    }
                }
            }
        }
    }

    return status;
}

} // namespace reduce
} // namespace coot