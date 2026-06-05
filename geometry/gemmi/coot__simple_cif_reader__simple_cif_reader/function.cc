#include "function.hh"
#include <gemmi/cif.hpp>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace coot {

void simple_cif_reader::simple_cif_reader_gemmi(const std::string &cif_dictionary_file_name) {
   names.clear();
   three_letter_codes.clear();

   std::filesystem::path p(cif_dictionary_file_name);
   if (!std::filesystem::exists(p) || !std::filesystem::is_regular_file(p)) {
      std::cout << "WARNIG:: cif dictionary " << cif_dictionary_file_name
                << " not found" << std::endl;
   } else {
      auto doc = gemmi::cif::read_file(cif_dictionary_file_name);
      if (doc.blocks.empty()) {
         std::cout << "Dirty mmCIF file? " << cif_dictionary_file_name
                   << std::endl;
         return;
      }
      for (auto & block : doc.blocks) {
         auto table = block.find_mmcif_category("_chem_comp");
         if (!table.ok()) continue;
         int n_col = table.find_column_position("name");
         int t_col = table.find_column_position("three_letter_code");
         if (n_col < 0 || t_col < 0) continue;
         for (auto row : table) {
            std::string name_val = row.at(n_col);
            std::string three_val = row.at(t_col);
            if (!name_val.empty() && !three_val.empty()) {
               names.push_back(name_val);
               three_letter_codes.push_back(three_val);
            }
         }
      }
   }
}

bool simple_cif_reader::has_restraints_for(const std::string &name) const {
   for (const auto &n : names) {
      if (n == name) return true;
   }
   return false;
}

} // namespace coot