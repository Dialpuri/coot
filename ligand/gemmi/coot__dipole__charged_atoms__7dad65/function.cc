#include "function.hh"

std::vector<std::pair<gemmi::Atom*, float>>
coot::dipole::charged_atoms_gemmi(
   std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::Residue*>> dict_res_pairs) {

   std::vector<std::pair<gemmi::Atom*, float>> charged_ats;

   for (size_t i = 0; i < dict_res_pairs.size(); i++) {
      gemmi::Residue* residue_p = dict_res_pairs[i].second;
      coot::dictionary_residue_restraints_t rest = dict_res_pairs[i].first;

      int n_dict_atom = (int)rest.atom_info.size();

      for (size_t i_res_at = 0; i_res_at < residue_p->atoms.size(); i_res_at++) {
         gemmi::Atom* at = &residue_p->atoms[i_res_at];
         std::string atom_name = at->name;
         for (int j = 0; j < n_dict_atom; j++) {
            if (rest.atom_info[j].partial_charge.first) {
               if (atom_name == rest.atom_info[j].atom_id_4c) {
                  std::pair<gemmi::Atom*, float> p(at, rest.atom_info[j].partial_charge.second);
                  charged_ats.push_back(p);
                  break;
               }
            }
         }
      }
   }
   return charged_ats;
}