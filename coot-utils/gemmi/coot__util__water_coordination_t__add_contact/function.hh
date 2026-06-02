#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>

// Ported contact_atom_t constructor (gemmi version)
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__contact_atoms_info_t__contact_atom_t__contact_atom_t__4d0713/gemmi/function.hh"

namespace coot {
namespace util {

// Gemmi-compatible contact info for a central atom (replaces coot::util::contact_atoms_info_t)
struct contact_info_gemmi {
    gemmi::Atom const* central;
    std::vector<coot::util::contact_atoms_info_t::contact_atom_t> contact_atoms;

    bool matches_atom(gemmi::Atom const* atom) const {
        return central == atom;
    }
};

// Gemmi-compatible water coordination (replaces coot::util::water_coordination_t)
// NOTE: This is a new gemmi-compatible struct, not the MMDB-based original.
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
        // The original condition: (alt_conf_1 == alt_conf_2) || (alt_conf_1 == "") || (alt_conf_2 == "")
        // Maps to: same altloc, or either has ' ' (no alt)
        if ((alt_conf_1 == alt_conf_2) || (alt_conf_1 == ' ') || (alt_conf_2 == ' ')) {
            // Filter out H water contacts.
            // In MMDB: atom->element == " H" (with leading space). In gemmi: element.name() == "H"
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
};

}  // namespace util
}  // namespace coot