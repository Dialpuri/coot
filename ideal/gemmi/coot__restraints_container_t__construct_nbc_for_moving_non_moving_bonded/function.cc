#include "function.hh"

namespace coot {

void construct_nbc_for_moving_non_moving_bonded_gemmi(
    std::vector<std::vector<unsigned int>>& filtered_non_bonded_atom_indices,
    unsigned int iat,
    unsigned int jat,
    const std::string& atom_name_1,
    const std::string& atom_name_2,
    const std::string& link_type,
    const coot::protein_geometry& geom)
{
  // Bounds check
  if (iat >= filtered_non_bonded_atom_indices.size()) return;

  dictionary_residue_link_restraints_t link = geom.link(link_type);

  if (!link.empty()) {
    bool add_it = true;

    for (unsigned int i = 0; i < link.link_bond_restraint.size(); i++) {
      if (atom_name_1 == link.link_bond_restraint[i].atom_id_1_4c() &&
          atom_name_2 == link.link_bond_restraint[i].atom_id_2_4c()) {
        add_it = false;
        break;
      }
    }
    for (unsigned int i = 0; i < link.link_angle_restraint.size(); i++) {
      if (atom_name_1 == link.link_angle_restraint[i].atom_id_1_4c() &&
          atom_name_2 == link.link_angle_restraint[i].atom_id_3_4c()) {
        add_it = false;
        break;
      }
    }
    for (unsigned int i = 0; i < link.link_torsion_restraint.size(); i++) {
      if (atom_name_1 == link.link_torsion_restraint[i].atom_id_1_4c() &&
          atom_name_2 == link.link_torsion_restraint[i].atom_id_4_4c()) {
        add_it = false;
        break;
      }
    }
    if (add_it) {
      filtered_non_bonded_atom_indices[iat].push_back(jat);
    }
  }
}

} // namespace coot