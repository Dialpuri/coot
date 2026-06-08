#include "function.hh"

namespace coot {

void restraints_container_t::set_z_occ_weights_gemmi() {
    atom_z_occ_weight.clear();
    atom_z_occ_weight.resize(n_atoms, 1.f);

    for (size_t i = 0; i < n_atoms; i++) {
        if (atom[i]) {
            float w = 1.f;

            // Occupancy weighting
            float occ = atom[i]->occ;
            if (occ < 1.f) {
                w = occ;
            }

            // Neutron refinement weighting (atomic number)
            if (do_neutron_refinement) {
                int z = atom[i]->element.atomic_number();
                if (z > 0) w = static_cast<float>(z);
            }

            // Cryo-EM weighting
            if (cryo_em_mode) {
                float b = atom[i]->b_iso;
                if (b > 0) {
                    float scale = 1.f / (1.f + b * 0.01f);
                    w *= scale;
                }
            }

            atom_z_occ_weight[i] = w;
        }
    }
}

void restraints_container_t::init_shared_post_gemmi(
    const std::vector<coot::atom_spec_t>& fixed_atom_specs)
{
  bonded_atom_indices.resize(n_atoms);

  set_has_hydrogen_atoms_state_gemmi();

  initial_position_params_vec.resize(3 * n_atoms);
  for (size_t i = 0; i < n_atoms; i++) {
    initial_position_params_vec[3 * i]     = atom[i]->pos.x;
    initial_position_params_vec[3 * i + 1] = atom[i]->pos.y;
    initial_position_params_vec[3 * i + 2] = atom[i]->pos.z;
  }

  // UDD is MMDB-specific; gemmi has no equivalent
  udd_bond_angle        = -1;
  udd_atom_index_handle = -1;

  use_map_gradient_for_atom.resize(n_atoms, false);
  if (!from_residue_vector) {
    for (size_t i = 0; i < n_atoms; i++) {
      gemmi::Residue* res = residue[i];
      if (res) {
        int seq = res->seqid.num.value;
        if (seq >= istart_res && seq <= iend_res) {
          if (!atom[i]->element.is_hydrogen())
            use_map_gradient_for_atom[i] = true;
        }
      }
    }
  } else {
    for (size_t i = 0; i < n_atoms; i++) {
      gemmi::Residue* res = residue[i];
      if (is_a_moving_residue_p_gemmi(res)) {
        if (!atom[i]->element.is_hydrogen() || do_hydrogen_atom_refinement)
          use_map_gradient_for_atom[i] = true;
      }
    }
  }

  set_z_occ_weights_gemmi();

  // assign_fixed_atom_indices — inline gemmi translation
  fixed_atom_indices.clear();
  for (const auto& spec : fixed_atom_specs) {
    for (size_t i = 0; i < n_atoms; i++) {
      if (spec.matches_spec(atom[i], residue[i], chain_ids[i])) {
        fixed_atom_indices.insert(static_cast<int>(i));
      }
    }
  }

  // Blank out fixed atoms from seeing electron density map gradients
  for (int idx : fixed_atom_indices)
    use_map_gradient_for_atom[static_cast<size_t>(idx)] = false;
}

} // namespace coot
