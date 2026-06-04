#include "function.hh"
#include <cstdlib>

namespace coot {

pepflip_using_difference_map::pepflip_using_difference_map(gemmi::Model* model)
    : model_(model)
{
}

std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth>>
pepflip_using_difference_map::make_random_other_pairs_gemmi(int n_others) const {

    std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth>> v;

    if (model_) {
        // Collect all atoms in file order (chain → residue → atom)
        std::vector<const gemmi::Atom*> atom_selection;
        for (const auto& chain : model_->chains) {
            for (const auto& res : chain.residues) {
                for (const auto& atom : res.atoms) {
                    atom_selection.push_back(&atom);
                }
            }
        }

        int n_selected_atoms = static_cast<int>(atom_selection.size());

        // Find all contacts: pairs (i,j) with min_dist <= distance <= max_dist
        // Mirrors mmdb::SeekContacts with min_dist=0.01, max_dist=3.0
        // SeekContacts returns both (i,j) and (j,i) for each contact pair
        std::vector<std::pair<int, int>> contacts;
        for (int i = 0; i < n_selected_atoms; i++) {
            for (int j = 0; j < n_selected_atoms; j++) {
                if (i == j) continue;
                double d = atom_selection[i]->pos.dist(atom_selection[j]->pos);
                if (d >= 0.01 && d <= 3.0) {
                    contacts.push_back({i, j});
                }
            }
        }

        int n_contacts = static_cast<int>(contacts.size());

        if (n_contacts > 0) {
            // Filter: select approximately n_others contacts
            // Original: filter = 2.0 * n_others / n_contacts
            // Then iterates over ALL contacts, checks ii<jj, calls rand() only for ii<jj
            float filter = 1.0f;
            if (n_contacts > n_others) {
                filter = 2.0f * static_cast<float>(n_others) / static_cast<float>(n_contacts);
            }
            float rmi = 1.0f / static_cast<float>(RAND_MAX);

            for (const auto& contact : contacts) {
                int ii = contact.first;
                int jj = contact.second;
                if (ii < jj) {
                    float f = static_cast<float>(std::rand()) * rmi;
                    if (f < filter) {
                        clipper::Coord_orth pt_1 = co_gemmi(atom_selection[ii]);
                        clipper::Coord_orth pt_2 = co_gemmi(atom_selection[jj]);
                        v.push_back({pt_1, pt_2});
                    }
                }
            }
        }
    }
    return v;
}

} // namespace coot