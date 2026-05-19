#pragma once
#include <gemmi/model.hpp>

namespace coot {

bool clashable_alt_confs_gemmi(const gemmi::Atom& at_1, const gemmi::Atom& at_2) {
    return gemmi::is_same_conformer(at_1.altloc, at_2.altloc);
}

} // namespace coot