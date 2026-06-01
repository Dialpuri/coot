#include "function.hh"
#include <iostream>

static std::string pad_atom_name_4(const std::string &name) {
    // MMDB stores atom names as 4-char PDB-style padded strings:
    //   left-padded with one space, right-padded with spaces to 4 chars.
    // e.g. "N" -> " N  ", "CA" -> " CA ", "CG1" -> "CG1 "
    std::string padded = ' ' + name;
    while (padded.size() < 4) {
        padded += ' ';
    }
    return padded;
}

std::vector<coot::atom_index_pair>
coot::monomer_utils::get_atom_index_pairs_gemmi(
  const std::vector<coot::atom_name_pair> &atom_name_pairs_in,
  const gemmi::Residue &residue) {

   int i_store_index;
   std::vector<coot::atom_index_pair> index_pairs;

   for (unsigned int ipair = 0; ipair < atom_name_pairs_in.size(); ipair++) {
      int ifound = 0;
      i_store_index = -1;
      for (int i = 0; i < static_cast<int>(residue.atoms.size()); i++) {
         std::string atomname = pad_atom_name_4(residue.atoms[i].name);
         if (atomname == atom_name_pairs_in[ipair].atom1) {
            i_store_index = i;
         }
      }
      if (i_store_index > -1) { // i.e. we found the first atom
         for (int i2 = 0; i2 < static_cast<int>(residue.atoms.size()); i2++) {
            std::string atomname = pad_atom_name_4(residue.atoms[i2].name);
            if (atomname == atom_name_pairs_in[ipair].atom2) {
               index_pairs.push_back(coot::atom_index_pair(i_store_index, i2));
            }
         }
      } else {
         std::cout << "first atom " << atom_name_pairs_in[ipair].atom1
                   << " not found in residue\n";
      }
   }

   if (index_pairs.size() != atom_name_pairs_in.size()) {
      std::cout << "failure to find all atom pair in residue atoms\n";
   } else {
      // std::cout << "DEBUG:: found all pairs in residue atoms\n";
   }
   return index_pairs;
}