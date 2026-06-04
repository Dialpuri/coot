#include "function.hh"
#include <string>
#include <vector>

namespace coot {
namespace util {

// Pad atom name to 4 characters the same way MMDB does
static std::string atom_name_to_4c(const std::string &name) {
  int ilen = static_cast<int>(name.length());
  if (ilen >= 4) return name.substr(0, 4);
  if (ilen == 1) return " " + name + "  ";
  if (ilen == 2) return " " + name + " ";
  if (ilen == 3) return " " + name;
  return name;
}

std::vector<std::vector<int>>
get_contact_indices_from_restraints_gemmi(gemmi::Residue *residue,
                                          const coot::dictionary_residue_restraints_t &restraints,
                                          bool regular_residue_flag,
                                          bool add_reverse_contacts) {

   int nResidueAtoms = static_cast<int>(residue->atoms.size());
   std::vector<std::vector<int>> contact_indices(nResidueAtoms);

   for (unsigned int ibr = 0; ibr < restraints.bond_restraint.size(); ibr++) {
      for (int iat = 0; iat < nResidueAtoms; iat++) {
         gemmi::Atom &atom_p = residue->atoms[iat];
         std::string at_name_4c(atom_name_to_4c(atom_p.name));
         if (restraints.bond_restraint[ibr].atom_id_1_4c() == at_name_4c) {
            int ibond_to = -1;  // initially unassigned
            std::string at_name_2;
            for (int iat2 = 0; iat2 < nResidueAtoms; iat2++) {
               gemmi::Atom &atom_p2 = residue->atoms[iat2];
               at_name_2 = atom_name_to_4c(atom_p2.name);
               if (restraints.bond_restraint[ibr].atom_id_2_4c() == at_name_2) {
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
   return contact_indices;
}

} // namespace util
} // namespace coot