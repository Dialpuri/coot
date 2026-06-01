#pragma once
#include <map>
#include <gemmi/model.hpp>
#include "geometry/hb-types.hh"

namespace coot {

class restraints_container_t_gemmi {
public:
    std::map<gemmi::Atom*, hb_t> H_atom_parent_energy_type_atom_map;

    bool H_parent_atom_is_donor_gemmi(gemmi::Atom* at) {
        bool state = false;
        std::map<gemmi::Atom*, hb_t>::const_iterator it;
        it = H_atom_parent_energy_type_atom_map.find(at);
        if (it != H_atom_parent_energy_type_atom_map.end()) {
            const hb_t& hbt = it->second;
            if (hbt == HB_DONOR || hbt == HB_BOTH)
                state = true;
        }
        return state;
    }
};

} // namespace coot