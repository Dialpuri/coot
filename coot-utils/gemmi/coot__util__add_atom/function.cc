#include "function.hh"
#include <gemmi/elem.hpp>
#include <iostream>

namespace coot {
namespace util {

namespace {
std::string trim_atom_name(const std::string& s) {
    size_t start = s.find_first_not_of(" \t");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t");
    return s.substr(start, end - start + 1);
}
} // anonymous namespace

bool add_atom_gemmi(
    gemmi::Residue* res,
    const std::string& atom_name_1,
    const std::string& atom_name_2,
    const std::string& atom_name_3,
    const std::string& alt_conf,
    double length,
    double angle,
    double torsion,
    const std::string& new_atom_name,
    const std::string& new_atom_ele,
    float new_atom_occ,
    float new_atom_b_factor
) {
    bool added_status = false;
    gemmi::Atom* a = nullptr;
    gemmi::Atom* b = nullptr;
    gemmi::Atom* c = nullptr;

    if (res) {
        std::string n1 = trim_atom_name(atom_name_1);
        std::string n2 = trim_atom_name(atom_name_2);
        std::string n3 = trim_atom_name(atom_name_3);

        for (gemmi::Atom& atom : res->atoms) {
            // alt_conf="" in original means no alternate location;
            // gemmi uses ' ' (space) for "no alt"
            bool alt_matches = false;
            if (alt_conf.empty()) {
                alt_matches = (atom.altloc == ' ' || atom.altloc == '\0');
            } else {
                alt_matches = (atom.altloc == alt_conf[0]);
            }

            if (alt_matches) {
                if (atom.name == n1) a = &atom;
                if (atom.name == n2) b = &atom;
                if (atom.name == n3) c = &atom;
            }
        }

        if (a && b && c) {
            clipper::Coord_orth ac(a->pos.x, a->pos.y, a->pos.z);
            clipper::Coord_orth bc(b->pos.x, b->pos.y, b->pos.z);
            clipper::Coord_orth cc(c->pos.x, c->pos.y, c->pos.z);
            double ang  = clipper::Util::d2rad(angle);
            double tors = clipper::Util::d2rad(torsion);
            clipper::Coord_orth pos(ac, bc, cc, length, ang, tors);

            gemmi::Atom new_atom;
            new_atom.name    = new_atom_name;
            new_atom.element = gemmi::Element(new_atom_ele);
            new_atom.pos     = gemmi::Position(pos.x(), pos.y(), pos.z());
            new_atom.altloc  = ' ';
            new_atom.occ     = new_atom_occ;
            new_atom.b_iso   = new_atom_b_factor;

            res->atoms.push_back(std::move(new_atom));
            added_status = true;
        } else {
            std::cout << "Failed to find all reference atoms : "
                      << atom_name_1 << " "
                      << atom_name_2 << " "
                      << atom_name_3 << std::endl;
        }
    }

    return added_status;
}

} // namespace util
} // namespace coot