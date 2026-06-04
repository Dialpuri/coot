#pragma once

#include <vector>
#include <gemmi/model.hpp>

namespace coot {
namespace util {
namespace map_to_model_correlation_stats_per_residue_run {

class residue_run_t_gemmi {
public:
   unsigned int idx_mid;
   unsigned int n_residues_per_blob;
   std::vector<const gemmi::Residue *> residues;

   residue_run_t_gemmi() {}

   explicit residue_run_t_gemmi(unsigned int n_residues_per_blob_in)
       : n_residues_per_blob(n_residues_per_blob_in) {
      idx_mid = n_residues_per_blob / 2;
   }

   explicit residue_run_t_gemmi(const std::vector<const gemmi::Residue *> &rr_in)
       : residues(rr_in) {
      n_residues_per_blob = rr_in.size();
      idx_mid = residues.size() / 2;
   }

   void add_residue(const gemmi::Residue *r) { residues.push_back(r); }

   const gemmi::Residue *residue_mid() const {
      if (residues.size() >= n_residues_per_blob) {
         return residues[idx_mid];
      } else {
         return nullptr;
      }
   }

   void add(const std::vector<const gemmi::Residue *> &rv) {
      residues.insert(residues.begin(), rv.cbegin(), rv.cend());
   }
};

} // namespace map_to_model_correlation_stats_per_residue_run
} // namespace util
} // namespace coot