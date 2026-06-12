#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <vector>
#include <string>
#include <iostream>

#include "mini-mol/mini-mol.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__minimol__atom__atom__9b7f96/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__minimol__residue__residue__7b6771/gemmi/function.hh"

namespace coot {
namespace minimol {

class molecule_gemmi {
public:
  std::vector<fragment> fragments;
  short have_cell;
  short have_spacegroup;

  molecule_gemmi(std::vector<gemmi::CRA> const& atom_selection,
                 int n_atoms,
                 std::vector<gemmi::Vec3> const& positions) {

    if (int(atom_selection.size()) != n_atoms) {
      std::cout << "ERROR:: inconsistence size in minimol molecule constructor"
                << std::endl;
      have_cell = 0;
      have_spacegroup = 0;
      return;
    }

    for (int iat = 0; iat < n_atoms; iat++) {

      gemmi::CRA cra = atom_selection[iat];
      gemmi::Atom* at = cra.atom;
      gemmi::Residue* residue_p = cra.residue;
      gemmi::Chain* chain_p = cra.chain;

      if (!at || !residue_p || !chain_p) {
        continue;
      }

      int resno = residue_p->seqid.num.value;
      std::string res_name = residue_p->name;
      std::string chain_id = chain_p->name;

      // now we have the properties of the atom, lets find where it
      // goes in the minimol.  First we need to find the chain, then
      // residue.  Then add the atom to the residue.

      bool found_fragment = false;
      int ifrag_for_atom = -1;
      for (unsigned int ifrag = 0; ifrag < fragments.size(); ifrag++) {
        if (fragments[ifrag].fragment_id == chain_id) {
          found_fragment = true;
          ifrag_for_atom = static_cast<int>(ifrag);
        }
      }

      if (!found_fragment) {
        fragment frag(chain_id);
        fragments.push_back(frag);
        ifrag_for_atom = static_cast<int>(fragments.size()) - 1;
      }

      // now find the residue, or add one.
      bool found_residue = false;
      if (resno <= fragments[ifrag_for_atom].max_residue_number()) {
        if (resno >= fragments[ifrag_for_atom].min_res_no()) {
          found_residue = true;
        }
      }

      if (!found_residue) {
        coot::minimol::residue res(resno);
        res.name = res_name;

        // Build a minimol::atom from gemmi atom data using the ported atom_gemmi
        coot::minimol::atom minimol_atom = coot::minimol::atom_gemmi(*at);
        minimol_atom.pos = clipper::Coord_orth(positions[iat].x, positions[iat].y, positions[iat].z);
        res.addatom(minimol_atom);
        try {
          fragments[ifrag_for_atom].addresidue(res, 1);
        }
        catch (const std::runtime_error &rte) {
          std::cout << "ERROR:: minimol constructor " << rte.what() << std::endl;
        }
      } else {
        coot::minimol::atom minimol_atom = coot::minimol::atom_gemmi(*at);
        minimol_atom.pos = clipper::Coord_orth(positions[iat].x, positions[iat].y, positions[iat].z);
        fragments[ifrag_for_atom][resno].addatom(minimol_atom);
      }
    }

    have_cell = 0;
    have_spacegroup = 0;
  }
};

}}  // namespace coot::minimol
