#include "function.hh"
#include <gemmi/contact.hpp>
#include <gemmi/neighbor.hpp>
#include <iostream>
#include <cmath>

namespace coot {
namespace high_res {

std::pair<clipper::Coord_orth, gemmi::Structure*>
get_middle_pos_gemmi(const gemmi::Structure& st) {
    std::pair<clipper::Coord_orth, gemmi::Structure*> r;

    if (st.models.empty()) {
        // Empty or null structure — return a non-null empty structure
        static gemmi::Structure empty_st;
        r.second = const_cast<gemmi::Structure*>(&empty_st);
        return r;
    }

    r.second = const_cast<gemmi::Structure*>(&st);

    gemmi::Model& model = const_cast<gemmi::Model&>(st.models[0]);
    double min_dist = 1.0;
    double max_dist = 15.0;

    // Collect all atoms with indices
    std::vector<gemmi::Atom*> atoms;
    for (gemmi::Chain& chain : model.chains)
        for (gemmi::Residue& res : chain.residues)
            for (gemmi::Atom& atom : res.atoms)
                atoms.push_back(&atom);

    int n_atoms = static_cast<int>(atoms.size());
    if (n_atoms == 0) return r;

    // Build neighbor search and find contacts
    gemmi::NeighborSearch ns(model, st.cell, max_dist);
    ns.populate(true);  // include hydrogens

    gemmi::ContactSearch cs(max_dist);
    cs.ignore = gemmi::ContactSearch::Ignore::Nothing;
    cs.twice = true;  // report both A-B and B-A, matching MMDB SeekContacts

    std::vector<gemmi::ContactSearch::Result> contacts = cs.find_contacts(ns);

    // Count contacts per atom (only partner1, matching original's id1 counting)
    std::vector<int> contact_count(n_atoms, 0);

    for (const auto& c : contacts) {
        double dist = std::sqrt(c.dist_sq);
        if (dist < min_dist) continue;
        if (!c.partner1.atom) continue;

        // Find index of partner1
        for (int i = 0; i < n_atoms; i++) {
            if (atoms[i] == c.partner1.atom) {
                contact_count[i]++;
                break;
            }
        }
    }

    // Find atom with most contacts
    int most_contacts_index = -1;
    int most_contacts = -1;
    for (int i = 0; i < n_atoms; i++) {
        if (contact_count[i] > most_contacts) {
            most_contacts_index = i;
            most_contacts = contact_count[i];
        }
    }

    if (most_contacts_index >= 0) {
        gemmi::Atom* at = atoms[most_contacts_index];
        r.first = clipper::Coord_orth(at->pos.x, at->pos.y, at->pos.z);
    }

    return r;
}

} // namespace high_res
} // namespace coot
