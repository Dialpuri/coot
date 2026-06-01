#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>

#include "geometry/protein-geometry.hh"

namespace coot {
namespace util {

namespace {

// Trim leading/trailing whitespace
inline std::string trim_ws(const std::string& s) {
  size_t a = s.find_first_not_of(' ');
  size_t b = s.find_last_not_of(' ');
  if (a == std::string::npos) return "";
  return s.substr(a, b - a + 1);
}

// Pad trimmed atom name to 4 chars (MMDB atom_id_4c style)
inline std::string atom_name_4c(const std::string& name) {
  std::string trimmed = trim_ws(name);
  int ilen = static_cast<int>(trimmed.length());
  if (ilen == 4) return trimmed;
  if (ilen == 1) return std::string(" ") + trimmed + std::string("  ");
  if (ilen == 2) return std::string(" ") + trimmed + std::string(" ");
  if (ilen == 3) return std::string(" ") + trimmed;
  return trimmed;
}

} // anonymous namespace

std::vector<std::vector<int>>
get_contact_indices_for_PRO_residue_gemmi(
    std::vector<const gemmi::Atom*> residue_atoms,
    coot::protein_geometry* geom_p) {

  int nResidueAtoms = static_cast<int>(residue_atoms.size());
  std::vector<std::vector<int>> contact_indices(nResidueAtoms);

  std::pair<bool, coot::dictionary_residue_restraints_t> rest =
    geom_p->get_monomer_restraints("PRO", coot::protein_geometry::IMOL_ENC_ANY);

  if (rest.first) {
    for (unsigned int ibr = 0; ibr < rest.second.bond_restraint.size(); ibr++) {
      for (int iat = 0; iat < nResidueAtoms; iat++) {
        const gemmi::Atom* atom_p = residue_atoms[iat];
        std::string at_name = atom_name_4c(atom_p->name);
        if (rest.second.bond_restraint[ibr].atom_id_1_4c() == at_name) {
          int ibond_to = -1;
          std::string at_name_2;
          for (int iat2 = 0; iat2 < nResidueAtoms; iat2++) {
            atom_p = residue_atoms[iat2];
            at_name_2 = atom_name_4c(atom_p->name);
            if (rest.second.bond_restraint[ibr].atom_id_2_4c() == at_name_2) {
              ibond_to = iat2;
              break;
            }
          }
          if (ibond_to != -1)
            contact_indices[iat].push_back(ibond_to);
        }
      }
    }
  }
  return contact_indices;
}

} // namespace util
} // namespace coot