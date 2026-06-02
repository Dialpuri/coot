#include "function.hh"
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__movable_atom/gemmi/function.hh"
#include <iostream>

namespace coot {
namespace molecule_t {

// Helper: find a matching atom in mol_model by spec from a CRA
static gemmi::Atom* find_atom_in_model(
    gemmi::Model& mol_model,
    const std::string& chain_id,
    int resno,
    char icode,
    const std::string& atom_name,
    char altloc)
{
    for (auto& ch : mol_model.chains) {
        if (ch.name != chain_id) continue;
        for (auto& res : ch.residues) {
            if (res.seqid.num.value != resno) continue;
            if (res.seqid.icode != icode) continue;
            for (auto& at : res.atoms) {
                if (at.name == atom_name && at.altloc == altloc) {
                    return &at;
                }
            }
        }
    }
    return nullptr;
}

void replace_coords_gemmi(
    gemmi::Model& mol_model,
    const std::vector<gemmi::CRA>& selection,
    bool change_altconf_occs_flag,
    bool replace_coords_with_zero_occ_flag)
{
    int n_atom = 0;

    for (const auto& cra : selection) {
        if (!cra.atom || !cra.residue || !cra.chain) continue;

        // Extract spec from the incoming CRA (the "new" atom)
        std::string chain_id = cra.chain->name;
        int resno = cra.residue->seqid.num.value;
        char icode = cra.residue->seqid.icode;
        std::string atom_name = cra.atom->name;
        char altloc = cra.atom->altloc;

        // Find matching atom in mol_model
        gemmi::Atom* mol_atom = find_atom_in_model(
            mol_model, chain_id, resno, icode, atom_name, altloc);

        if (!mol_atom) {
            std::cout << "WARNING:: bad atom idx -1 (atom not found in model)" << std::endl;
            continue;
        }

        if (change_altconf_occs_flag) {
            // "change alt conf occs" mode — use the new atom's occupancy
            if (coot::movable_atom_gemmi(mol_atom, replace_coords_with_zero_occ_flag)) {
                mol_atom->pos = cra.atom->pos;
                mol_atom->occ = cra.atom->occ;
                // keep mol_atom->b_iso unchanged
                n_atom++;
            }
        } else {
            // "don't change alt confs" mode — keep original occupancy
            if (coot::movable_atom_gemmi(mol_atom, replace_coords_with_zero_occ_flag)) {
                mol_atom->pos = cra.atom->pos;
                // keep mol_atom->occ and mol_atom->b_iso unchanged
                n_atom++;
            }
        }
    }
}

} // namespace molecule_t
} // namespace coot