#include "function.hh"

namespace coot {
namespace molecule_t {

void add_target_position_restraint_gemmi(
    const std::string &atom_cid,
    float pos_x, float pos_y, float pos_z,
    gemmi::Structure &st,
    std::vector<std::pair<gemmi::Atom *, clipper::Coord_orth>> &atoms_with_position_restraints)
{
    gemmi::Atom *at = cid_to_atom_gemmi(atom_cid, st);
    if (at) {
        bool done = false;
        for (auto &entry : atoms_with_position_restraints) {
            if (entry.first == at) {
                entry.second = clipper::Coord_orth(pos_x, pos_y, pos_z);
                done = true;
                break;
            }
        }
        if (!done) {
            clipper::Coord_orth p(pos_x, pos_y, pos_z);
            atoms_with_position_restraints.push_back(std::make_pair(at, p));
        }
    }
}

} // namespace molecule_t
} // namespace coot