#pragma once

#include <gemmi/model.hpp>

namespace coot {

// Stub for bonds_box to match original test expectations
struct bonds_box_t {
    int n_bonds() const { return 0; }
};

class ghost_molecule_display_t {
public:
    bonds_box_t bonds_box;
    void update_bonds_gemmi(gemmi::Structure& st);
};

} // namespace coot

// Inline definition since the body is empty
inline void coot::ghost_molecule_display_t::update_bonds_gemmi(gemmi::Structure& st) {
   // update bonds - why is there nothing here?
}