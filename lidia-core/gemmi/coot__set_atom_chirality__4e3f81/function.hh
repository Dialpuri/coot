#pragma once

#include <gemmi/model.hpp>
#include <rdkit/GraphMol/RDKitBase.h>
#include "geometry/protein-geometry.hh"
#include "utils/coot-utils.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__get_chiral_tag/gemmi/function.hh"
#include <string>

namespace coot {

inline void set_atom_chirality_gemmi(
    RDKit::Atom* rdkit_at,
    const gemmi::Atom* atom_p,
    const gemmi::Residue* residue_p,
    const dictionary_residue_restraints_t& restraints) {

  // set the chirality
  // (if this atom has restraints-style chiral info)
  bool done_chiral = false;

  // Normalize: gemmi's atom.name is not padded to 4 chars like MMDB
  std::string atom_name_stripped = coot::util::remove_whitespace(atom_p->name);

  for (unsigned int ichi = 0; ichi < restraints.chiral_restraint.size(); ichi++) {
    const dict_chiral_restraint_t& cr = restraints.chiral_restraint[ichi];
    std::string cr_atom_c_stripped = coot::util::remove_whitespace(cr.atom_id_c_4c());
    if (cr_atom_c_stripped == atom_name_stripped) {
      done_chiral = true;
      if (!cr.has_unassigned_chiral_volume()) {
        rdkit_at->setProp("mmcif_chiral_N1", coot::util::remove_whitespace(cr.atom_id_1_4c()));
        rdkit_at->setProp("mmcif_chiral_N2", coot::util::remove_whitespace(cr.atom_id_2_4c()));
        rdkit_at->setProp("mmcif_chiral_N3", coot::util::remove_whitespace(cr.atom_id_3_4c()));
        if (!cr.is_a_both_restraint()) {
          // e.g. RDKit::Atom::CHI_TETRAHEDRAL_CCW;
          RDKit::Atom::ChiralType chiral_tag =
              get_chiral_tag_gemmi(residue_p, restraints, atom_p);
          rdkit_at->setChiralTag(chiral_tag);

          std::string bc("positive");
          if (cr.volume_sign == dict_chiral_restraint_t::CHIRAL_RESTRAINT_NEGATIVE)
            bc = "negative";

          rdkit_at->setProp("mmcif_chiral_volume_sign", bc);

        } else {
          std::string bc("both");
          rdkit_at->setProp("mmcif_chiral_volume_sign", bc);
        }
      }
    }
  }

  // set chirality
  // (if this atom has Chemical Component Dictionary style chirality (R/S pdbx_stereo_config_flag)
  //
  if (!done_chiral) {
    for (unsigned int i = 0; i < restraints.atom_info.size(); i++) {
      std::string atom_info_name_stripped = coot::util::remove_whitespace(restraints.atom_info[i].atom_id_4c);
      if (atom_info_name_stripped == atom_name_stripped) {
        set_atom_chirality(rdkit_at, restraints.atom_info[i]);
      }
    }
  }
}

} // namespace coot
