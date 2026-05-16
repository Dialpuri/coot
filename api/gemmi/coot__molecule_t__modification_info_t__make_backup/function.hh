#pragma once

#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <gemmi/pdb.hpp>
#include <gemmi/to_pdb.hpp>

namespace coot {

struct modification_info_t {
   struct save_info_t {
      std::string file_name;
      std::string modification_info_string;
      save_info_t() = default;
      save_info_t(const std::string &fn, const std::string &mi)
         : file_name(fn), modification_info_string(mi) {}
   };
   std::vector<save_info_t> save_info;
   int modification_index = 0;
   std::string backup_dir = "coot-backup";
   std::string mol_name;
   bool is_mmcif_flag = false;

   std::string get_backup_file_name_from_index(int index) const {
      return backup_dir + "/test-" + std::to_string(index) + ".pdb";
   }

   // Port of make_backup from MMDB to gemmi.
   // Writes a PDB backup of `st` to disk, appends to save_info,
   // and sets modification_index.  Returns "null molecule" if st
   // has no models or no atoms — in that case info is unchanged.
   std::string make_backup_gemmi(gemmi::Structure &st,
                                 const std::string &modification_info_string)
   {
      if (st.models.empty())
         return "null molecule";
      size_t n = 0;
      for (auto &m : st.models)
         for (auto &c : m.chains)
            for (auto &r : c.residues)
               n += r.atoms.size();
      if (n == 0)
         return "null molecule";

      int idx = static_cast<int>(save_info.size());
      std::string fn = get_backup_file_name_from_index(idx);

      std::filesystem::create_directories(backup_dir);

      std::ofstream ofs(fn);
      gemmi::write_pdb(st, ofs);
      ofs.close();

      save_info.emplace_back(fn, modification_info_string);
      modification_index = static_cast<int>(save_info.size());
      return "";
   }
};

} // namespace coot