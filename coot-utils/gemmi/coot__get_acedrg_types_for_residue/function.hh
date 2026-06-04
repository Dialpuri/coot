#pragma once

#include <gemmi/model.hpp>
#include "coot-utils/acedrg-types-for-residue.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

acedrg_types_for_residue_t
get_acedrg_types_for_residue_gemmi(const gemmi::Residue& residue, int imol_enc,
                                   const protein_geometry& geom) {
  acedrg_types_for_residue_t types;
  std::string residue_type = residue.name;
  auto r = geom.get_monomer_restraints(residue_type, imol_enc);
  if (r.first) {
    const auto& restraints = r.second;
    for (unsigned int ib = 0; ib < restraints.bond_restraint.size(); ib++) {
      const auto& bond_restraint = restraints.bond_restraint[ib];
      const std::string& atom_name_1 = bond_restraint.atom_id_1();
      const std::string& atom_name_2 = bond_restraint.atom_id_2();

      int idx_1 = -1;
      int idx_2 = -1;

      // Build a local vector of atom pointers from the residue
      std::vector<const gemmi::Atom*> residue_atoms;
      for (auto& atom : residue.atoms) {
        residue_atoms.push_back(&atom);
      }

      for (int iat = 0; iat < static_cast<int>(residue_atoms.size()); iat++) {
        const gemmi::Atom* at = residue_atoms[iat];
        // gemmi stores unpadded names; MMDB restraints use 4-char padded names
        std::string atom_name = pad_atom_name(at->name, at->element.name());
        if (atom_name == atom_name_1) idx_1 = iat;
        if (atom_name == atom_name_2) idx_2 = iat;
      }

      if (idx_1 != -1) {
        if (idx_2 != -1) {
          clipper::Coord_orth at_pos_1 = coot::co_gemmi(residue_atoms[idx_1]);
          clipper::Coord_orth at_pos_2 = coot::co_gemmi(residue_atoms[idx_2]);
          double bb = (at_pos_2 - at_pos_1).lengthsq();
          double bond_length = std::sqrt(bb);

          std::string type_1;
          std::string type_2;

          for (unsigned int ii = 0; ii < restraints.atom_info.size(); ii++) {
            const auto& atom = restraints.atom_info[ii];
            if (atom.atom_id_4c == atom_name_1) type_1 = atom.acedrg_atom_type;
            if (atom.atom_id_4c == atom_name_2) type_2 = atom.acedrg_atom_type;
          }

          if (!type_1.empty()) {
            if (!type_2.empty()) {
              bool in_same_ring_flag = restraints.in_same_ring(atom_name_1, atom_name_2);
              acedrg_types_for_bond_t bt(atom_name_1, atom_name_2, type_1, type_2, bond_length, in_same_ring_flag);
              types.bond_types.push_back(bt);
            }
          }
        }
      }
    }
  }
  return types;
}

} // namespace coot