#include "function.hh"
#include <gemmi/neighbor.hpp>

namespace coot {

void atom_overlaps_container_t::fill_ligand_atom_neighbour_map_gemmi() {
    if (!model || !res_central) {
        return;
    }

    double max_dist = 2.3;
    gemmi::UnitCell uc; // default (non-periodic)

    gemmi::NeighborSearch ns(*model, uc, max_dist);
    ns.populate(true); // include hydrogens

    // Collect atoms from res_central
    std::vector<gemmi::Atom *> residue_atoms;
    for (gemmi::Atom &atom : res_central->atoms) {
        residue_atoms.push_back(&atom);
    }

    for (int i = 0; i < static_cast<int>(residue_atoms.size()); i++) {
        gemmi::Atom *at = residue_atoms[i];
        auto neighbors = ns.find_atoms(at->pos, '\0', 0.0, max_dist);

        for (auto *mark : neighbors) {
            gemmi::CRA cra = mark->to_cra(*model);
            if (cra.atom && cra.residue == res_central && cra.atom != at) {
                double radius = get_vdw_radius_ligand_atom_gemmi(cra.atom);
                ligand_atom_neighbour_map[i].push_back({cra.atom, radius});
            }
        }
    }
}

double atom_overlaps_container_t::get_vdw_radius_ligand_atom_gemmi(gemmi::Atom *at) {
    if (at && at->element.name()[0] != '\0' && at->element.vdw_r() > 0) {
        return at->element.vdw_r();
    }
    return 1.7;
}

} // namespace coot