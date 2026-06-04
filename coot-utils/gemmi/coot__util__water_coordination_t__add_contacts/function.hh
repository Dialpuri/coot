#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <vector>
#include <string>
#include <cmath>

// Ported contact_atom_t constructor (gemmi version)
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__contact_atoms_info_t__contact_atom_t__contact_atom_t__4d0713/gemmi/function.hh"

namespace coot {
namespace util {

// Gemmi-compatible contact info for a central atom
struct contact_info_gemmi {
    gemmi::Atom const* central;
    std::vector<coot::util::contact_atoms_info_t::contact_atom_t> contact_atoms;

    bool matches_atom(gemmi::Atom const* atom) const {
        return central == atom;
    }
};

// Gemmi-compatible water coordination
class water_coordination_t {
public:
    std::vector<contact_info_gemmi> atom_contacts;

    void add_contact_gemmi(
        gemmi::Atom const* atom_central,
        gemmi::Atom const* atom_contactor,
        double const mat[4][4])
    {
        char alt_conf_1 = atom_contactor->altloc;
        char alt_conf_2 = atom_central->altloc;

        // In MMDB, empty altLoc ("") means no alternate location.
        // In gemmi, altloc ' ' (space) means no alternate location.
        if ((alt_conf_1 == alt_conf_2) || (alt_conf_1 == ' ') || (alt_conf_2 == ' ')) {
            // Filter out H water contacts.
            if (atom_contactor->element.name() != "H") {
                contact_atoms_info_t::contact_atom_t con_at =
                    contact_atoms_info_t::contact_atom_t::contact_atom_t_gemmi(
                        atom_contactor, atom_central, mat);

                bool found = false;
                for (size_t i = 0; i < atom_contacts.size(); ++i) {
                    if (atom_contacts[i].matches_atom(atom_central)) {
                        atom_contacts[i].contact_atoms.push_back(con_at);
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    contact_info_gemmi cai;
                    cai.central = atom_central;
                    cai.contact_atoms.push_back(con_at);
                    atom_contacts.push_back(cai);
                }
            }
        }
    }

    void add_contacts_gemmi(
        std::vector<gemmi::Atom const*> const& water_selection,
        std::vector<gemmi::Atom const*> const& atom_selection,
        double min_dist, double max_dist,
        double const mat[4][4])
    {
        double min_dist_sq = min_dist * min_dist;
        double max_dist_sq = max_dist * max_dist;

        for (auto const* water_atom : water_selection) {
            for (auto const* selected_atom : atom_selection) {
                gemmi::Vec3 diff = water_atom->pos - selected_atom->pos;
                double dist_sq = diff.length_sq();
                if (dist_sq >= min_dist_sq && dist_sq <= max_dist_sq) {
                    add_contact_gemmi(water_atom, selected_atom, mat);
                }
            }
        }
    }
};

}  // namespace util
}  // namespace coot