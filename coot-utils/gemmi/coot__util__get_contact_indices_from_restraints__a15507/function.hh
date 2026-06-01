#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "geometry/protein-geometry.hh"
#include <vector>
#include <string>

namespace coot {
namespace util {

// Pad atom name to 4 characters the same way atom_id_mmdb_expand does
static std::string atom_name_to_4c(const std::string &name) {
  int ilen = static_cast<int>(name.length());
  if (ilen == 4) return name;
  if (ilen == 1) {
    return " " + name + "  ";
  }
  if (ilen == 2) {
    return " " + name + " ";
  }
  if (ilen == 3) {
    return " " + name;
  }
  return name;
}

std::vector<std::vector<int>>
get_contact_indices_from_restraints_gemmi(gemmi::Residue *residue,
                                          const coot::protein_geometry *geom_p,
                                          bool regular_residue_flag,
                                          bool add_reverse_contacts) {

  int nResidueAtoms = static_cast<int>(residue->atoms.size());
  std::vector<std::vector<int>> contact_indices(nResidueAtoms);
  std::string restype(residue->name);

  int n_monomers = static_cast<int>(geom_p->size());

  auto residue_has_deuterium_atoms = [](gemmi::Residue *residue) {
    int nResidueAtoms = static_cast<int>(residue->atoms.size());
    bool has_deuterium_atoms = false;
    for (int iat = 0; iat < nResidueAtoms; iat++) {
      gemmi::Atom &atom_p = residue->atoms[iat];
      std::string atom_ele(atom_p.element.name());
      if (atom_ele == "D") {
        has_deuterium_atoms = true;
        break;
      }
    }
    return has_deuterium_atoms;
  };

  bool has_deuterium_atoms = residue_has_deuterium_atoms(residue);

  for (int icomp = 0; icomp < n_monomers; icomp++) {
    const dictionary_residue_restraints_t &dict = (*geom_p)[icomp].second;
    if (dict.residue_info.comp_id == restype) {
      for (unsigned int ibr = 0; ibr < dict.bond_restraint.size(); ibr++) {
        for (int iat = 0; iat < nResidueAtoms; iat++) {
          gemmi::Atom &atom_p = residue->atoms[iat];
          std::string at_name_4c(atom_name_to_4c(atom_p.name));

          if (dict.bond_restraint[ibr].atom_id_1_4c() == at_name_4c) {
            int ibond_to = -1;
            std::string at_name_2_4c;
            for (int iat2 = 0; iat2 < nResidueAtoms; iat2++) {
              gemmi::Atom &atom_p2 = residue->atoms[iat2];
              at_name_2_4c = atom_name_to_4c(atom_p2.name);
              if (dict.bond_restraint[ibr].atom_id_2_4c() == at_name_2_4c) {
                ibond_to = iat2;
                break;
              }
            }
            if (ibond_to > -1) {
              if (add_reverse_contacts == 0) {
                if (regular_residue_flag) {
                  contact_indices[iat].push_back(ibond_to);
                } else {
                  contact_indices[ibond_to].push_back(iat);
                }
              } else {
                contact_indices[ibond_to].push_back(iat);
                contact_indices[iat].push_back(ibond_to);
              }
            }
          }

          if (has_deuterium_atoms) {
            std::string dict_atom_name_1 = dict.bond_restraint[ibr].atom_id_1_4c();
            if (dict_atom_name_1[0] == 'H') dict_atom_name_1[0] = 'D';
            if (dict_atom_name_1.size() > 1 && dict_atom_name_1[1] == 'H') dict_atom_name_1[1] = 'D';
            if (dict_atom_name_1 == at_name_4c) {
              int ibond_to = -1;
              for (int iat2 = 0; iat2 < nResidueAtoms; iat2++) {
                gemmi::Atom &atom_p2 = residue->atoms[iat2];
                std::string at_name_2_4c = atom_name_to_4c(atom_p2.name);
                std::string dict_atom_name_2 = dict.bond_restraint[ibr].atom_id_2_4c();
                if (dict_atom_name_2[0] == 'H') dict_atom_name_2[0] = 'D';
                if (dict_atom_name_2.size() > 1 && dict_atom_name_2[1] == 'H') dict_atom_name_2[1] = 'D';
                if (dict_atom_name_2 == at_name_2_4c) {
                  ibond_to = iat2;
                  break;
                }
              }
              if (ibond_to > -1) {
                if (add_reverse_contacts == 0) {
                  if (regular_residue_flag) {
                    contact_indices[iat].push_back(ibond_to);
                  } else {
                    contact_indices[ibond_to].push_back(iat);
                  }
                } else {
                  contact_indices[ibond_to].push_back(iat);
                  contact_indices[iat].push_back(ibond_to);
                }
              }
            }
          }
        }
      }
    }
  }
  return contact_indices;
}

} // namespace util
} // namespace coot