#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <gemmi/cifdoc.hpp>

namespace coot {
namespace protein_geometry {

struct dict_link_plane_restraint_t {
  std::string plane_id;
  std::vector<std::string> atom_ids;
  std::vector<int> atom_comp_ids;
  std::vector<double> dist_esds;

  std::string atom_id(size_t n) const { return atom_ids[n]; }
  double dist_esd() const { return dist_esds[0]; }
  std::size_t n_atoms() const { return atom_ids.size(); }
};

struct dictionary_residue_link_restraints_t {
  std::string link_id;
  std::vector<dict_link_plane_restraint_t> link_plane_restraint;
};

inline void link_add_plane(std::string const &link_id,
                           std::string const &atom_id,
                           std::string const &plane_id,
                           int atom_comp_id,
                           double dist_esd,
                           std::vector<dictionary_residue_link_restraints_t> &dict_link_res_restraints) {
  dict_link_plane_restraint_t new_plane_restraint;
  new_plane_restraint.plane_id = plane_id;
  new_plane_restraint.dist_esds.push_back(dist_esd);
  new_plane_restraint.atom_comp_ids.push_back(atom_comp_id);
  new_plane_restraint.atom_ids.push_back(atom_id);

  bool found = false;
  for (auto &res_entry : dict_link_res_restraints) {
    if (res_entry.link_id == link_id) {
      res_entry.link_plane_restraint.push_back(new_plane_restraint);
      found = true;
      break;
    }
  }
  if (!found) {
    dictionary_residue_link_restraints_t new_entry;
    new_entry.link_id = link_id;
    new_entry.link_plane_restraint.push_back(new_plane_restraint);
    dict_link_res_restraints.push_back(new_entry);
  }
}

inline void link_plane_gemmi(gemmi::cif::Loop &loop,
                             std::vector<dictionary_residue_link_restraints_t> &dict_link_res_restraints) {

   for (size_t j = 0; j < loop.length(); j++) {

      std::string link_id;
      std::string atom_id, plane_id;
      double dist_esd = 0.0;
      int atom_comp_id = 0;

      int ierr_tot = 0;

      int col_link_id = loop.find_tag("link_id");
      if (col_link_id >= 0)
         link_id = gemmi::cif::as_string(loop.val(j, col_link_id));
      else
         ierr_tot++;

      int col_atom_id = loop.find_tag("atom_id");
      if (col_atom_id >= 0)
         atom_id = gemmi::cif::as_string(loop.val(j, col_atom_id));
      else
         ierr_tot++;

      int col_atom_comp = loop.find_tag("atom_comp_id");
      if (col_atom_comp >= 0)
         atom_comp_id = gemmi::cif::as_int(loop.val(j, col_atom_comp));
      else
         ierr_tot++;

      int col_plane_id = loop.find_tag("plane_id");
      if (col_plane_id >= 0)
         plane_id = gemmi::cif::as_string(loop.val(j, col_plane_id));
      else
         ierr_tot++;

      int col_dist_esd = loop.find_tag("dist_esd");
      if (col_dist_esd >= 0)
         dist_esd = std::stod(loop.val(j, col_dist_esd));
      else
         ierr_tot++;

      if (ierr_tot == 0) {
         link_add_plane(link_id, atom_id, plane_id, atom_comp_id, dist_esd, dict_link_res_restraints);
      } else {
         std::cout << "problem reading link plane mmCIFLoop" << std::endl;
      }
   }
}

// Simulates what coot::protein_geometry::init_standard() populates for link restraints
inline std::vector<dictionary_residue_link_restraints_t> init_standard_link_restraints() {
  std::vector<dictionary_residue_link_restraints_t> result;
  for (int i = 0; i < 108; i++) {
    dictionary_residue_link_restraints_t entry;
    entry.link_id = "STANDARD_LINK_" + std::to_string(i);
    result.push_back(entry);
  }
  return result;
}

} // namespace protein_geometry
} // namespace coot