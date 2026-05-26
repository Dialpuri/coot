#pragma once
#include <vector>
#include <utility>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <clipper/clipper.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__extents__0761f9/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__centre_of_molecule/gemmi/function.hh"

namespace coot {
namespace util {

inline std::vector<std::pair<double, double>>
spherically_averaged_molecule_gemmi(const gemmi::Structure& st, float angstroms_per_bin, const std::string& chain_id = "") {
    std::vector<std::pair<double, double>> vp;

    std::pair<clipper::Coord_orth, clipper::Coord_orth> e = coot::util::extents_gemmi(st, chain_id);
    std::pair<bool, clipper::Coord_orth> cc = coot::centre_of_molecule_gemmi(st);
    if (!cc.first) return vp;
    clipper::Coord_orth c = cc.second;

    double diag_len_sqrd = (e.second - e.first).lengthsq();
    double diag_len = std::sqrt(diag_len_sqrd);
    double radius_max = 0.5 * diag_len;

    int n_bins = static_cast<int>(radius_max / angstroms_per_bin) + 1;
    vp.resize(n_bins);
    for (int ibin = 0; ibin < n_bins; ibin++) {
        vp[ibin].first = (static_cast<float>(ibin) + 0.5) * angstroms_per_bin;
        vp[ibin].second = 0.0;
    }

    // Traverse atoms
    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            if (!chain_id.empty() && chain.name != chain_id) continue;
            for (const gemmi::Residue& res : chain.residues) {
                for (const gemmi::Atom& at : res.atoms) {
                    clipper::Coord_orth co(at.pos.x, at.pos.y, at.pos.z);
                    float dist = std::sqrt((co - c).lengthsq());
                    int bin_id = static_cast<int>(dist / angstroms_per_bin);
                    if (bin_id >= n_bins) {
                        std::cout << "ERROR:: bin error! " << std::endl;
                    } else {
                        vp[bin_id].second += 1.0;
                    }
                }
            }
        }
    }

    return vp;
}

} // namespace util
} // namespace coot