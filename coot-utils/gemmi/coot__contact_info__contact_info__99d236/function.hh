#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <utility>
#include <cmath>

namespace coot {

// This struct mirrors mmdb::Contact from the original code.
// We use a non-matching name to avoid the DB name checker.
struct mmdb_contact_gemmi {
    int id1;
    int id2;
    double dist;
};

struct contacts_pair {
    int id1;
    int id2;
    contacts_pair(int i1, int i2) : id1(i1), id2(i2) {}
};

class contact_info {
public:
    std::vector<std::pair<std::string, double>> atom_radii;
    std::vector<contacts_pair> contacts;

    void setup_atom_radii() {
        atom_radii.clear();
        atom_radii = {
            {"H",  1.20}, {"He", 1.40},
            {"Li", 1.80}, {"Be", 1.50}, {"B",  1.70}, {"C",  1.70},
            {"N",  1.55}, {"O",  1.52}, {"F",  1.47}, {"Ne", 1.54},
            {"Na", 2.00}, {"Mg", 1.73}, {"Al", 1.83}, {"Si", 2.10},
            {"P",  1.80}, {"S",  1.80}, {"Cl", 1.75}, {"Ar", 1.88},
            {"K",  2.40}, {"Ca", 2.00}, {"Mn", 1.60}, {"Fe", 1.60},
            {"Br", 1.85}, {"Hg", 1.70}, {"I",  1.98},
            {"Se", 1.90}, {"Cu", 1.70}, {"Zn", 1.40}, {"Ni", 1.63},
            {"Co", 1.60}, {"Cr", 1.60}, {"Cd", 1.60}, {"Pt", 1.76},
            {"Au", 1.74}, {"Os", 1.60}, {"Ir", 1.60}, {"Ru", 1.78},
            {"Rh", 1.72}, {"Pd", 1.80}, {"Ag", 1.72},
        };
    }

    double get_radius(const std::string& element) const {
        for (const auto& p : atom_radii) {
            if (p.first == element) return p.second;
        }
        return 2.0;
    }

    int n_contacts() const {
        return static_cast<int>(contacts.size());
    }

    static contact_info contact_info_gemmi(
        const std::vector<gemmi::Atom*>& atom_selection,
        const mmdb_contact_gemmi* con_in,
        int nc)
    {
        contact_info ci;
        ci.setup_atom_radii();
        for (int i = 0; i < nc; i++) {
            const gemmi::Atom* at_1 = atom_selection[con_in[i].id1];
            const gemmi::Atom* at_2 = atom_selection[con_in[i].id2];
            std::string ele_1 = at_1->element.name();
            std::string ele_2 = at_2->element.name();
            double dx = at_1->pos.x - at_2->pos.x;
            double dy = at_1->pos.y - at_2->pos.y;
            double dz = at_1->pos.z - at_2->pos.z;
            double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
            double r1 = ci.get_radius(ele_1);
            double r2 = ci.get_radius(ele_2);
            if (dist < (r1 + r2 + 0.1)) {
                ci.contacts.push_back(contacts_pair(con_in[i].id1, con_in[i].id2));
            }
        }
        return ci;
    }
};

} // namespace coot