#include "function.hh"
#include "coords/gemmi/closest_approach/function.hh"
#include <gemmi/model.hpp>
#include <iostream>

// Gemmi port of coot::high_res::fill_globular_protein
void coot_high_res_fill_globular_protein_gemmi(
    coot::minimol::molecule &globular_mol,
    const coot::minimol::molecule &mol,
    const clipper::Coord_orth &target_pos_in,
    const gemmi::Structure &structure)
{
    clipper::Coord_orth sum_atoms = target_pos_in;
    clipper::Coord_orth target_pos = target_pos_in;

    // Lets add 20 atoms at the target_pos;
    sum_atoms = clipper::Coord_orth(20.0 * target_pos.x(),
                                    20.0 * target_pos.y(),
                                    20.0 * target_pos.z());
    double n_atoms = 20;

    clipper::Coord_orth t;

    globular_mol.set_cell_symm(mol);

    for (unsigned int ifrag = 0; ifrag < mol.fragments.size(); ifrag++) {
        int igfrag = globular_mol.fragment_for_chain(mol[ifrag].fragment_id);

        coot::minimol::residue residue(1, "ALA");

        for (int ires = mol[ifrag].min_res_no();
             ires <= mol[ifrag].max_residue_number();
             ires++)
        {
            for (unsigned int iat = 0;
                 iat < mol[ifrag][ires].n_atoms();
                 iat++)
            {
                // Use gemmi version of closest_approach
                t = closest_approach_gemmi(
                    mol[ifrag][ires][iat].pos,
                    target_pos,
                    structure);

                sum_atoms += t;
                n_atoms += 1.0;

                residue.addatom(" C  ", " C", t, "", 1.0, 30.0);
            }
        }

        try {
            globular_mol[igfrag].addresidue(residue, 0);
        }
        catch (const std::runtime_error &rte) {
            std::cout << "ERROR:: fill_globular_protein() "
                      << rte.what() << std::endl;
        }
    }

    std::cout << "DEBUG:: ##################### globular_molecule created"
              << std::endl;
}