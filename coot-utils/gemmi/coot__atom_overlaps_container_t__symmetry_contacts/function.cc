#include "function.hh"

namespace coot {

std::vector<atom_overlap_t>
atom_overlaps_container_t::symmetry_contacts_gemmi(float d) {

  std::vector<atom_overlap_t> v;

  if (!mol_ || !res_central_) {
    return v;
  }

  // Get the first model
  if (mol_->models.empty()) {
    return v;
  }
  gemmi::Model& model = mol_->models[0];

  // Get unit cell and space group
  gemmi::UnitCell cell = mol_->cell;
  const gemmi::SpaceGroup* sg_ptr = mol_->find_spacegroup();
  if (!sg_ptr) {
    return v;
  }
  gemmi::GroupOps ops = sg_ptr->operations();
  int n_symm = ops.order();

  // Find the chain that contains res_central_
  gemmi::Chain* chain_of_res = nullptr;
  for (gemmi::Chain& ch : model.chains) {
    for (gemmi::Residue& res : ch.residues) {
      if (&res == res_central_) {
        chain_of_res = &ch;
        break;
      }
    }
    if (chain_of_res) break;
  }

  if (!chain_of_res) {
    return v;
  }

  // Build selection of central residue atoms using the existing gemmi port
  coot::residue_spec_t spec(chain_of_res->name, res_central_->seqid.num.value);
  std::vector<gemmi::Atom*> atom_selection_res = coot::select_atoms_gemmi(spec, mol_);

  if (atom_selection_res.empty()) {
    return v;
  }

  // Collect ALL atoms in the model
  std::vector<gemmi::Atom*> atom_selection_all;
  for (gemmi::Chain& ch : model.chains) {
    for (gemmi::Residue& res : ch.residues) {
      for (gemmi::Atom& atom : res.atoms) {
        atom_selection_all.push_back(&atom);
      }
    }
  }

  int shift_lim = 2;
  double min_contact_dist = static_cast<double>(d);

  // For each combination of symmetry op and translation shift
  // skip the identity (i=0 with shifts 0,0,0)
  for (int x_shift = -shift_lim; x_shift <= shift_lim; x_shift++) {
    for (int y_shift = -shift_lim; y_shift <= shift_lim; y_shift++) {
      for (int z_shift = -shift_lim; z_shift <= shift_lim; z_shift++) {
        for (int i_symm = 0; i_symm < n_symm; i_symm++) {
          if (!(x_shift == 0 && y_shift == 0 && z_shift == 0 && i_symm == 0)) {
            gemmi::Op symop = ops.get_op(i_symm);

            // For each atom in the full model, compute its symmetry-equivalent position
            for (gemmi::Atom* at_all : atom_selection_all) {
              gemmi::Position orig_pos = at_all->pos;

              // Convert to fractional coordinates
              gemmi::Fractional frac = cell.fractionalize(orig_pos);

              // Apply the rotation part of the symmetry operation to fractional coords
              double nx = static_cast<double>(symop.rot[0][0]) * frac.x
                        + static_cast<double>(symop.rot[0][1]) * frac.y
                        + static_cast<double>(symop.rot[0][2]) * frac.z
                        + static_cast<double>(symop.tran[0])
                        + static_cast<double>(x_shift);

              double ny = static_cast<double>(symop.rot[1][0]) * frac.x
                        + static_cast<double>(symop.rot[1][1]) * frac.y
                        + static_cast<double>(symop.rot[1][2]) * frac.z
                        + static_cast<double>(symop.tran[1])
                        + static_cast<double>(y_shift);

              double nz = static_cast<double>(symop.rot[2][0]) * frac.x
                        + static_cast<double>(symop.rot[2][1]) * frac.y
                        + static_cast<double>(symop.rot[2][2]) * frac.z
                        + static_cast<double>(symop.tran[2])
                        + static_cast<double>(z_shift);

              // Convert back to orthogonal (Cartesian) coordinates
              gemmi::Position new_pos = cell.orthogonalize(gemmi::Fractional(nx, ny, nz));

              // Check distance to each central residue atom
              for (gemmi::Atom* at_res : atom_selection_res) {
                double dist_sq = (at_res->pos - new_pos).length_sq();
                if (dist_sq < min_contact_dist * min_contact_dist) {
                  v.emplace_back(at_res, at_all);
                }
              }
            }
          }
        }
      }
    }
  }

  return v;
}

} // namespace coot
