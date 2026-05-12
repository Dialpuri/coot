#pragma once

#include <gemmi/elem.hpp>
#include <gemmi/model.hpp>

namespace coot {

struct atom_overlaps_container_t {
    bool is_ss_bonded_or_CYS_CYS_SGs_gemmi(const gemmi::Atom& at_1, const gemmi::Residue& res_1,
                                           const gemmi::Atom& at_2, const gemmi::Residue& res_2) const {
        if (at_1.element != gemmi::Element("S") || at_2.element != gemmi::Element("S"))
            return false;
        if (res_1.name != "CYS" || res_2.name != "CYS")
            return false;
        if (at_1.name != "SG" || at_2.name != "SG")
            return false;
        return true;
    }
};

} // namespace coot