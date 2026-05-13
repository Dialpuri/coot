#include "function.hh"

namespace coot {

float molecule_t_get_temperature_factor_of_atom_gemmi(const std::string &cid, gemmi::Structure &st) {
    float b = -1.1f;
    gemmi::Atom *at = molecule_t::cid_to_atom_gemmi(cid, st);
    if (at) {
        b = at->b_iso;
    }
    return b;
}

} // namespace coot