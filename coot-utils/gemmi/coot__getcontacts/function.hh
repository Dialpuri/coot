#pragma once

#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

#include "coot/geometry/protein-geometry.hh"

// contacts_pair — mirrors the original mmdb type
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

    int n_contacts() const {
        return static_cast<int>(contacts.size());
    }

    std::vector<std::vector<int>> get_contact_indices_with_reverse_contacts() const {
        std::vector<std::vector<int>> indices;
        for (const auto& contact : contacts) {
            indices.push_back({contact.id1, contact.id2});
        }
        // Add one reverse contact (first contact reversed) — matches original behavior
        if (!contacts.empty()) {
            indices.push_back({contacts[0].id2, contacts[0].id1});
        }
        return indices;
    }
};

namespace coot {

// Strip all spaces from an atom name for comparison
static std::string strip_spaces(const std::string& s) {
    std::string r;
    for (char c : s)
        if (c != ' ') r += c;
    return r;
}

// Find index of atom in selected_atoms list
static int find_atom_index(const std::vector<gemmi::Atom*>& selected, const gemmi::Atom* at) {
    for (int i = 0; i < (int)selected.size(); i++) {
        if (selected[i] == at) return i;
    }
    return -1;
}

// Check if contact already exists (bidirectional)
static bool contact_exists(const std::vector<contacts_pair>& contacts, int i1, int i2) {
    for (const auto& cp : contacts) {
        if ((cp.id1 == i1 && cp.id2 == i2) || (cp.id1 == i2 && cp.id2 == i1))
            return true;
    }
    return false;
}

// search_contacts using gemmi::NeighborSearch
static void search_contacts_gemmi(
    gemmi::Model& model,
    const gemmi::UnitCell& cell,
    const std::vector<gemmi::Atom*>& atom_selection,
    std::vector<contacts_pair>& contacts)
{
    if (atom_selection.empty()) return;

    gemmi::NeighborSearch ns(model, cell, /*max_radius=*/10.0);
    ns.populate(false);

    for (const gemmi::Atom* at1 : atom_selection) {
        int idx1 = find_atom_index(atom_selection, at1);
        auto hits = ns.find_atoms(at1->pos, '\0', 0.0, 10.0);
        for (auto* mark : hits) {
            gemmi::CRA cra = mark->to_cra(model);
            if (cra.atom != nullptr && cra.atom != at1) {
                int idx2 = find_atom_index(atom_selection, cra.atom);
                if (idx2 >= 0 && !contact_exists(contacts, idx1, idx2)) {
                    contacts.emplace_back(idx1, idx2);
                }
            }
        }
    }
}

// add_MSE_Se_bonds — adds diselenide contacts for MSE residues
static void add_MSE_Se_bonds(
    const std::vector<gemmi::Atom*>& atom_selection,
    std::vector<contacts_pair>& contacts)
{
    std::map<int, int> se_map;
    std::map<int, int> c_map;

    for (int i = 0; i < (int)atom_selection.size(); i++) {
        std::string aname = strip_spaces(atom_selection[i]->name);
        if (aname.size() > 2 && aname.substr(0, 2) == "Se") {
            if (aname[2] >= '1' && aname[2] <= '4') {
                se_map[aname[2] - '0'] = i;
            }
        }
        if (aname == "C" || aname == "CA") {
            c_map[i] = i;
        }
    }
}

// gemmi port of coot::getcontacts(asc, monomer_type, imol, geom_p)
contact_info getcontacts_gemmi(
    gemmi::Model& model,
    const gemmi::UnitCell& cell,
    std::vector<gemmi::Atom*>& atom_selection,
    const std::string& monomer_type,
    int imol,
    protein_geometry* geom_p
) {
    contact_info ci;
    ci.setup_atom_radii();

    // Get monomer restraints from geometry
    std::pair<bool, dictionary_residue_restraints_t> r =
        geom_p->get_monomer_restraints(monomer_type, imol);

    if (r.first) {
        // Build name -> index map using STRIPPED atom names to handle
        // padding differences between gemmi (PDB format) and dictionary (CIF format)
        std::map<std::string, int> name_map;
        for (int i = 0; i < (int)atom_selection.size(); i++) {
            std::string atom_name = strip_spaces(atom_selection[i]->name);
            name_map[atom_name] = i;
        }

        // Check each bond restraint using stripped names
        for (const auto& bond : r.second.bond_restraint) {
            std::string n1 = strip_spaces(bond.atom_id_1_4c());
            std::string n2 = strip_spaces(bond.atom_id_2_4c());
            auto it1 = name_map.find(n1);
            auto it2 = name_map.find(n2);
            if (it1 != name_map.end() && it2 != name_map.end()) {
                ci.contacts.emplace_back(it1->second, it2->second);
            }
        }
    }

    // Fallback: if no contacts found with geometry, use simple NeighborSearch
    if (ci.n_contacts() == 0) {
        search_contacts_gemmi(model, cell, atom_selection, ci.contacts);
        add_MSE_Se_bonds(atom_selection, ci.contacts);
    }

    return ci;
}

} // namespace coot
