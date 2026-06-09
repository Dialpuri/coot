#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <clipper/clipper.h>
#include "mini-mol/mini-mol.hh"

namespace coot { namespace minimol {

inline molecule molecule_gemmi(
    const std::vector<gemmi::CRA>& atom_selection,
    int n_atoms,
    const std::vector<gemmi::Position>& atoms)
{
    molecule mol;

    if (int(atoms.size()) != n_atoms) {
        std::cout << "ERROR:: inconsistence size in minimol molecule constructor"
                  << std::endl;
        return mol;
    }

    for (int iat = 0; iat < n_atoms; iat++) {

        gemmi::CRA const& cra = atom_selection[iat];
        gemmi::Chain const* chain_p = cra.chain;
        gemmi::Residue const* residue_p = cra.residue;
        gemmi::Atom const* atom_p = cra.atom;

        int resno = residue_p->seqid.num.value;
        std::string res_name = residue_p->name;
        std::string chain_id = chain_p->name;

        // Find the fragment (chain) or create one
        bool found_fragment = false;
        int ifrag_for_atom = -1;
        for (unsigned int ifrag = 0; ifrag < mol.fragments.size(); ifrag++) {
            if (mol.fragments[ifrag].fragment_id == chain_id) {
                found_fragment = true;
                ifrag_for_atom = static_cast<int>(ifrag);
            }
        }

        if (!found_fragment) {
            fragment frag(chain_id);
            mol.fragments.push_back(frag);
            ifrag_for_atom = static_cast<int>(mol.fragments.size()) - 1;
        }

        // Find the residue or add one
        bool found_residue = false;
        if (resno <= mol.fragments[ifrag_for_atom].max_residue_number()) {
            if (resno >= mol.fragments[ifrag_for_atom].min_res_no()) {
                found_residue = true;
            }
        }

        // Construct minimol atom from gemmi atom
        coot::minimol::atom minimol_atom;
        if (atom_p) {
            minimol_atom.name = atom_p->name;
            minimol_atom.element = atom_p->element.name();
            minimol_atom.altLoc = atom_p->altloc;
            minimol_atom.occupancy = atom_p->occ;
            minimol_atom.temperature_factor = atom_p->b_iso;
            minimol_atom.int_user_data = -1;
        }
        // Override position from the atoms vector (as in original)
        minimol_atom.pos = clipper::Coord_orth(atoms[iat].x, atoms[iat].y, atoms[iat].z);

        if (!found_residue) {
            coot::minimol::residue res(resno);
            res.name = res_name;
            res.addatom(minimol_atom);
            try {
                mol.fragments[ifrag_for_atom].addresidue(res, 1);
            } catch (const std::runtime_error &rte) {
                std::cout << "ERROR:: minimol constructor " << rte.what() << std::endl;
            }
        } else {
            mol.fragments[ifrag_for_atom][resno].addatom(minimol_atom);
        }
    }

    return mol;
}

}} // namespace coot::minimol
