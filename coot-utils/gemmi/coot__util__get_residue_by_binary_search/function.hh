#pragma once
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot { namespace util {

inline gemmi::Residue *
get_residue_by_binary_search_gemmi(const std::string &chain_id,
                                   int res_no, const std::string &insertion_code,
                                   gemmi::Structure &st) {
   gemmi::Residue *res = nullptr;
   bool found_res = false;

   // Convert insertion code: MMDB uses "" for no insertion code, gemmi uses ' '
   char ins_code_char = insertion_code.empty() ? ' ' : insertion_code[0];

   for (gemmi::Model& model : st.models) {
      for (gemmi::Chain& chain : model.chains) {
         if (chain.name != chain_id) continue;

         // Binary search setup
         int nres = static_cast<int>(chain.residues.size());
         if (nres == 0) continue;

         int top_idx = nres - 1;
         int bottom_idx = 0;
         int idx_trial_prev = -999;

         while (!found_res) {
            int idx_delta = top_idx - bottom_idx;
            int idx_trial = bottom_idx + idx_delta / 2;

            if (idx_trial < 0 || idx_trial >= nres) break;

            gemmi::Residue &residue_this = chain.residues[idx_trial];

            if (residue_this.seqid.num.value == res_no) {
               char this_ins_code = residue_this.seqid.icode;
               if (this_ins_code == ins_code_char) {
                  res = &residue_this;
                  found_res = true;
                  break;
               }
            }

            if (!found_res) {
               if (top_idx == bottom_idx) break;
               if (idx_trial == bottom_idx) break;

               if (residue_this.seqid.num.value > res_no) {
                  top_idx = idx_trial;
               }
               if (residue_this.seqid.num.value < res_no) {
                  bottom_idx = idx_trial;
               }
               if (idx_trial == idx_trial_prev) break;
            }

            idx_trial_prev = idx_trial;
         }

         if (!found_res) {
            // Fallback: linear search
            for (int ires = 0; ires < nres; ++ires) {
               gemmi::Residue &residue_p = chain.residues[ires];
               if (residue_p.seqid.num.value == res_no) {
                  if (residue_p.seqid.icode == ins_code_char) {
                     res = &residue_p;
                     found_res = true;
                     break;
                  }
               }
            }
         }

         if (found_res) break;
      }
      if (found_res) break;
   }

   return res;
}

}}