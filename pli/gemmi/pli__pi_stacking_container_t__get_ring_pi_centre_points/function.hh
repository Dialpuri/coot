#pragma once

#include <string>
#include <vector>
#include <utility>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>

namespace pli {

class pi_stacking_container_t {
public:
    // Original MMDB version:
    // std::pair<clipper::Coord_orth, clipper::Coord_orth>
    // get_ring_pi_centre_points(const std::vector<std::string> &ring_atom_names,
    //                           mmdb::Residue *res_ref) const;

    // Gemmi version - uses gemmi::Residue instead of mmdb::Residue
    // Note: gemmi::Residue doesn't have parent pointer, so we need to pass chain name
    std::pair<clipper::Coord_orth, clipper::Coord_orth>
    get_ring_pi_centre_points_gemmi(const std::vector<std::string> &ring_atom_names,
                                    const gemmi::Residue &res_ref,
                                    const std::string &chain_name) const;

private:
    // Helper function - can throw an exception if not enough points found in pts.
    std::pair<clipper::Coord_orth, clipper::Coord_orth>
    ring_centre_and_normal(const std::vector<clipper::Coord_orth> &pts) const;
};

} // namespace pli