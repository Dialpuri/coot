#include "function.hh"

namespace pli {

double find_water_protein_length_gemmi(const gemmi::Residue* ligand_residue,
                                       const gemmi::Model* model) {
    double dist = 100.0;
    double dist_sqrd = dist * dist;
    double dist_sqrd_init = dist_sqrd;

    // Get ligand residue atoms
    const auto& ligand_atoms = ligand_residue->atoms;

    for (const gemmi::Chain& chain : model->chains) {
        for (const gemmi::Residue& residue : chain.residues) {
            if (&residue != ligand_residue) {
                // Skip water molecules
                if (residue.name == "HOH")
                    continue;

                // Check each atom pair
                for (const gemmi::Atom& lig_atom : ligand_atoms) {
                    for (const gemmi::Atom& res_atom : residue.atoms) {
                        std::string ele = res_atom.element.name();
                        if (ele == "O" || ele == "N") {
                            clipper::Coord_orth pt_1(lig_atom.pos.x,
                                                     lig_atom.pos.y,
                                                     lig_atom.pos.z);
                            clipper::Coord_orth pt_2(res_atom.pos.x,
                                                     res_atom.pos.y,
                                                     res_atom.pos.z);
                            double d2 = (pt_1 - pt_2).lengthsq();
                            if (d2 < dist_sqrd) {
                                dist_sqrd = d2;
                            }
                        }
                    }
                }
            }
        }
    }

    if (dist_sqrd < dist_sqrd_init) {
        dist = std::sqrt(dist_sqrd);
    }

    return dist;
}

} // namespace pli