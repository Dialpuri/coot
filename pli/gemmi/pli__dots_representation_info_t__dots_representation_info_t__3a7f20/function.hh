#pragma once

#include "pli/dots-representation-info.hh"
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace pli {

// gemmi port of the constructor dots_representation_info_t(mmdb::Manager*, mmdb::Manager*)
// The original only sets is_closed=0 and does commented-out selection work;
// the default constructor already sets is_closed=false, so we return a default
// constructed object.
inline dots_representation_info_t dots_representation_info_t_gemmi(
    gemmi::Structure const* mol,
    gemmi::Structure const* mol_exclude) {
    dots_representation_info_t result;
    (void)mol;
    (void)mol_exclude;
    return result;
}

} // namespace pli