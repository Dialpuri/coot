#include "function.hh"
#include <algorithm>
#include "coot-utils/coot-coord-utils.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

std::pair<bool, double> coot::angle_betwen_plane_and_vector_gemmi(
    gemmi::Residue *residue_p,
    const std::vector<std::string> &ring_atom_names,
    const std::string &altconf_in,
    const clipper::Coord_orth &vector) {

    std::pair<bool, double> r(false, 0);

    if (residue_p) {
        std::vector<clipper::Coord_orth> ring_atom_positions;
        for (gemmi::Atom &at : residue_p->atoms) {
            std::string atom_name = at.name;
            std::string alt_conf(1, at.altloc);
            std::vector<std::string>::const_iterator it =
                std::find(ring_atom_names.begin(), ring_atom_names.end(), atom_name);
            if (it != ring_atom_names.end()) {
                if (alt_conf == altconf_in) {
                    clipper::Coord_orth pos = coot::co_gemmi(&at);
                    ring_atom_positions.push_back(pos);
                }
            }

            if (ring_atom_positions.size() > 4) {
                coot::lsq_plane_info_t l(ring_atom_positions);
                double angle = l.angle(vector);
                r.first = true;
                r.second = angle;
            }
        }
    }

    return r;
}