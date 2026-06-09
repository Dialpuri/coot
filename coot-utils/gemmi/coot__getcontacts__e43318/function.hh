#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

#include "geometry/protein-geometry.hh"

// Include the verified gemmi port of contact_info constructor (with alt_conf)
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__contact_info__contact_info__e43318/gemmi/function.hh"

namespace coot {

/// gemmi port of coot::getcontacts(asc, alt_conf, monomer_type, imol, geom_p)
inline coot::contact_info::contact_info_result_gemmi getcontacts_gemmi(
    const std::vector<gemmi::Atom*>& atom_selection,
    const std::string& alt_conf,
    const std::string& monomer_type,
    int imol,
    coot::protein_geometry* geom_p)
{
    coot::contact_info::contact_info_result_gemmi ci = coot::contact_info::contact_info_gemmi(
        atom_selection, alt_conf, monomer_type, imol, geom_p);

    if (ci.n_contacts() == 0) {
        // Fallback: return empty result (equivalent to coot::getcontacts(asc)
        // which is not yet ported, and only matters for empty selections)
        return coot::contact_info::contact_info_result_gemmi{};
    }
    return ci;
}

} // namespace coot
