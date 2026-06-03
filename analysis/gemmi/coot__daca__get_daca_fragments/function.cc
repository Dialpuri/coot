#include "function.hh"
#include "analysis/daca.hh"
#include <algorithm>

namespace {

std::string trim_ws(const std::string& s) {
    auto start = s.find_first_not_of(' ');
    if (start == std::string::npos) return std::string();
    auto end = s.find_last_not_of(' ');
    return s.substr(start, end - start + 1);
}

} // anonymous namespace

std::vector<std::vector<const gemmi::Atom*>>
coot::get_daca_fragments_gemmi(const gemmi::Residue* reference_residue_p) {

   std::vector<std::vector<const gemmi::Atom*>> v;
   std::string res_name(reference_residue_p->name);
   coot::daca daca_obj;
   std::vector<std::vector<std::string>> atom_name_vec_vec =
       daca_obj.atom_names_for_fragments(res_name);

   for (const auto& atom_names : atom_name_vec_vec) {
      std::vector<const gemmi::Atom*> atom_vec;
      for (const auto& atom_name : atom_names) {
         std::string trimmed_target = trim_ws(atom_name);
         for (const auto& atom : reference_residue_p->atoms) {
            std::string trimmed_atom = trim_ws(atom.name);
            if (trimmed_target == trimmed_atom) {
               if (atom.altloc == '\0') {
                  atom_vec.push_back(&atom);
               }
            }
         }
      }
      if (atom_names.size() == atom_vec.size()) {
         v.push_back(atom_vec);
      }
   }
   return v;
}