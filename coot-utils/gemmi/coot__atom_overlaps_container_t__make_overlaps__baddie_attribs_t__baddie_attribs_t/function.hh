#pragma once
#include <gemmi/model.hpp>
#include <string>
#include <vector>

namespace coot {

struct atom_overlaps_container_t {
    struct make_overlaps {
        struct baddie_attribs_t {
            gemmi::Atom* cr_at;
            gemmi::Atom* n_at;
            float r_1;
            float r_2;
            float d;
            float o;
            bool is_hydrogen_bond;
            bool hydrogen_atom_is_first_atom;
            
            baddie_attribs_t(gemmi::Atom* cr_at, gemmi::Atom* n_at, float r_1, float r_2, float d, float o,
                             bool is_hydrogen_bond, bool hydrogen_atom_is_first_atom)
                : cr_at(cr_at), n_at(n_at), r_1(r_1), r_2(r_2), d(d), o(o),
                  is_hydrogen_bond(is_hydrogen_bond), hydrogen_atom_is_first_atom(hydrogen_atom_is_first_atom) {}
            
            static bool sorter(const baddie_attribs_t &b1, const baddie_attribs_t &b2) {
                return b2.o < b1.o;
            }
        };
    };
};

} // namespace coot