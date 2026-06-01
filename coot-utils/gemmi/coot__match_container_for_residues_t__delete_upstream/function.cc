#include "function.hh"
#include <algorithm>
#include <set>

namespace coot {

void
match_container_for_residues_t::delete_upstream_gemmi(gemmi::Structure& st, bool from_first,
                                                       const std::vector<gemmi::CRA>& atoms0,
                                                       const std::vector<gemmi::CRA>& atoms1) {

    // Pick the selection based on from_first
    const std::vector<gemmi::CRA>& sel_atoms = from_first ? atoms0 : atoms1;

    if (sel_atoms.empty())
        return;

    // All atoms in the selection come from the first model
    gemmi::Model& model = st.models[0];

    // Find the chain from the first selected atom
    const gemmi::CRA& first_cra = sel_atoms[0];
    if (first_cra.chain == nullptr)
        return;

    gemmi::Chain& chain_ref = *first_cra.chain;

    // Collect residues to delete (use a set of pointers for dedup)
    std::set<gemmi::Residue*> delete_these_residues;

    // Find the "matcher" residue — the residue containing the first atom in sel_atoms
    // that is part of atom_pairs
    gemmi::Residue* matchers_residue = nullptr;
    bool found_matchers = false;

    for (const gemmi::CRA& cra : sel_atoms) {
        gemmi::Atom* at = cra.atom;
        if (at == nullptr) continue;

        for (unsigned int ip = 0; ip < atom_pairs.size(); ip++) {
            if (from_first) {
                if (atom_pairs[ip].first.atom == at) {
                    found_matchers = true;
                    matchers_residue = atom_pairs[ip].first.residue;
                    break;
                }
            } else {
                if (atom_pairs[ip].second.atom == at) {
                    found_matchers = true;
                    matchers_residue = atom_pairs[ip].second.residue;
                    break;
                }
            }
        }

        if (found_matchers)
            break;

        // If not found, and this atom's residue != matchers_residue,
        // add it to delete_these_residues (if not already present)
        if (cra.residue != matchers_residue && cra.residue != nullptr) {
            delete_these_residues.insert(cra.residue);
        }
    }

    // Now delete the collected residues
    if (!delete_these_residues.empty()) {
        // Use std::remove_if to remove all at once (avoid pointer invalidation)
        auto it = std::remove_if(chain_ref.residues.begin(),
                                 chain_ref.residues.end(),
                                 [&delete_these_residues, matchers_residue](gemmi::Residue& r) {
                                     gemmi::Residue* rp = &r;
                                     return delete_these_residues.count(rp) && (rp != matchers_residue);
                                 });
        chain_ref.residues.erase(it, chain_ref.residues.end());
    }
}

} // namespace coot