#include "function.hh"

#include "coot-utils/coot-coord-utils.hh"

namespace coot {
namespace util {

missing_atom_info_gemmi missing_atoms_gemmi(
    gemmi::Structure* st,
    bool do_missing_hydrogen_atoms_flag,
    coot::protein_geometry* geom_p
) {

   bool ignore_missing_OXT = true;
   bool ignore_missing_OP3 = true;

   std::vector<gemmi::Residue*> residues_with_missing_atoms;
   std::vector<std::string> residues_no_dictionary;
   std::map<gemmi::Residue*, std::vector<std::string>> residue_missing_atom_names_map;
   std::vector<std::pair<gemmi::Residue*, std::vector<gemmi::Atom*>>> atoms_in_coords_but_not_in_dict;

   if (st) {

      std::vector<coot::util::dict_residue_atom_info_t> residue_atoms;

      // Use model 0 (gemmi is 0-indexed)
      gemmi::Model& model = (*st).models[0];

      for (gemmi::Chain& chain : model.chains) {
         for (gemmi::Residue& residue : chain.residues) {
            std::string residue_name(residue.name);

            short int found_dict = 0;
            std::vector<coot::util::dict_atom_info_t> residue_dict_atoms;

            for (unsigned int idict = 0; idict < residue_atoms.size(); idict++) {
               std::string tmp_str = residue_atoms[idict].residue_name;
               if (residue_name == tmp_str) {
                  residue_dict_atoms = residue_atoms[idict].atom_info;
                  found_dict = 1;
               }
            }

            if (!found_dict) {
               coot::util::dict_residue_atom_info_t residue_atoms_for_a_residue(residue_name, geom_p);
               if (!residue_atoms_for_a_residue.is_empty_p()) {
                  residue_atoms.push_back(residue_atoms_for_a_residue);
                  residue_dict_atoms = residue_atoms_for_a_residue.atom_info;
                  found_dict = 1;
               }
            }

            if (!found_dict) {
               residues_no_dictionary.push_back(residue_name);
            } else {

               std::vector<coot::util::dict_atom_info_t> dict_atom_names_pairs;
               for (unsigned int iat = 0; iat < residue_dict_atoms.size(); iat++) {
                  if ((do_missing_hydrogen_atoms_flag == false) &&
                      (residue_dict_atoms[iat].is_Hydrogen_flag == 1)) {
                     // do nothing — skip hydrogens when not requested
                  } else {
                     coot::util::dict_atom_info_t p(residue_dict_atoms[iat].name, 0);
                     bool really_missing = true;

                     if (ignore_missing_OXT) {
                        if (residue_dict_atoms[iat].name == " OXT")
                           really_missing = false;
                     }
                     if (ignore_missing_OP3) {
                        if (residue_dict_atoms[iat].name == " OP3")
                           really_missing = false;
                     }

                     if (really_missing) {
                        dict_atom_names_pairs.push_back(p);
                     }
                  }
               }

               // For every atom in the residue, check if it's in the dictionary
               for (gemmi::Atom& atom : residue.atoms) {
                  // gemmi atom.name is unpadded (e.g. "CA"), but dictionary names
                  // are padded to 4 chars (e.g. " CA "). Use pad_atom_name to match.
                  std::string atom_name_padded = coot::pad_atom_name(atom.name, atom.element.name());
                  for (unsigned int idictat = 0; idictat < dict_atom_names_pairs.size(); idictat++) {
                     if (atom_name_padded == dict_atom_names_pairs[idictat].name) {
                        dict_atom_names_pairs[idictat].is_Hydrogen_flag = 1; // mark as found
                        break;
                     }
                  }
               }

               std::vector<std::string> missing_atom_names;
               for (unsigned int idictat = 0; idictat < dict_atom_names_pairs.size(); idictat++) {
                  if (!dict_atom_names_pairs[idictat].is_Hydrogen_flag) {
                     missing_atom_names.push_back(dict_atom_names_pairs[idictat].name);
                  }
               }

               if (!missing_atom_names.empty()) {
                  residues_with_missing_atoms.push_back(&residue);
                  residue_missing_atom_names_map[&residue] = missing_atom_names;
               }
            }
         }
      }
   }

   missing_atom_info_gemmi mai(residues_no_dictionary,
                               residues_with_missing_atoms,
                               atoms_in_coords_but_not_in_dict);
   mai.residue_missing_atom_names_map = residue_missing_atom_names_map;

   return mai;
}

} // namespace util
} // namespace coot
