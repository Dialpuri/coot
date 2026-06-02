#pragma once

#include <vector>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__lidia_utils__co/gemmi/function.hh"

namespace coot {
namespace chem_feat_clust {

struct ligand_t {
    gemmi::Residue* residue;
};

struct chem_feat_clust_t {
    std::vector<ligand_t> ligands;
};

inline std::vector<clipper::Coord_orth> get_ligands_coords_gemmi(const chem_feat_clust_t& clust) {
    std::vector<clipper::Coord_orth> v;
    for (const auto& lig : clust.ligands) {
        if (lig.residue) {
            for (const auto& atom : lig.residue->atoms) {
                v.push_back(coot::lidia_utils::co_gemmi(&atom));
            }
        }
    }
    return v;
}

} // namespace chem_feat_clust
} // namespace coot