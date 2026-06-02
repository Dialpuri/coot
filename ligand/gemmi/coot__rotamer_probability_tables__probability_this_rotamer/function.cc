#include "ligand/rotamer.hh"
#include <vector>
#include <string>
#include <utility>
#include <stdexcept>

namespace coot {

rotamer_probability_info_t probability_this_rotamer_gemmi(
  const rotamer_probability_tables &tables_obj,
  unsigned int i_table,
  const std::vector<std::pair<int,float>> &chi_angles)
{
   const a_rotamer_table &table = tables_obj[i_table];

   // Inline chi_angles_to_bins logic (no _gemmi port exists)
   std::vector<int> bins(chi_angles.size());
   for (size_t i = 0; i < chi_angles.size(); i++) {
      int chi_number = chi_angles[i].first;
      float chi_angle = chi_angles[i].second;

      // Normalize angle to [-180, 180]
      while (chi_angle > 180.0f) chi_angle -= 360.0f;
      while (chi_angle < -180.0f) chi_angle += 360.0f;

      int n_samples = 0;
      if (chi_number == 1) n_samples = table.n_chi1_samples_per_360;
      else if (chi_number == 2) n_samples = table.n_chi2_samples_per_360;
      else if (chi_number == 3) n_samples = table.n_chi3_samples_per_360;
      else if (chi_number == 4) n_samples = table.n_chi4_samples_per_360;

      if (n_samples == 0) {
         std::string mess = "ERROR: chi number " + std::to_string(chi_number) +
             " is not valid for " + table.residue_name;
         throw std::runtime_error(mess);
      }

      float bin_size = 360.0f / static_cast<float>(n_samples);
      int bin = static_cast<int>(chi_angle / bin_size);
      if (bin < 0) bin += n_samples;
      bins[i] = bin;
   }

   if (bins.size() != chi_angles.size()) {
      throw std::runtime_error("ERROR:: bin size and chi_angles size do not match");
   }

   unsigned int n = table.n_chis;
   if (bins.size() != n) {
      std::string mess = "ERROR:: not enough chi angles found. bin.size() (";
      mess += std::to_string(bins.size());
      mess += ") and n chis do not match ";
      mess += " (should be ";
      mess += std::to_string(n);
      mess += ") for ";
      mess += table.residue_name;
      throw std::runtime_error(mess);
   }

   float pr = 0.0f;
   if (n == 1) pr = table.pr_chi_1[bins[0]];
   if (n == 2) pr = table.pr_chi_1_2[bins[0]][bins[1]];
   if (n == 3) pr = table.pr_chi_1_2_3[bins[0]][bins[1]][bins[2]];
   if (n == 4) pr = table.pr_chi_1_2_3_4[bins[0]][bins[1]][bins[2]][bins[3]];

   if (n < 1 || n > 4) {
      std::string mess = "ERROR: bad nbins chis " + std::to_string(n);
      mess += " for i_table ";
      mess += std::to_string(i_table);
      throw std::runtime_error(mess);
   }

   return rotamer_probability_info_t(rotamer_probability_info_t::OK, pr * 100.0f, table.residue_name);
}

} // namespace coot