#include <cmath>
#include "function.hh"

static int round_up_by_hundreds(int num) {
   float a = static_cast<float>(num + 99) * 0.01f;
   float f = floorf(a);
   int ii = static_cast<int>(f) * 100;
   return ii;
}

namespace coot {
namespace cho {

std::pair<short int, int>
next_residue_number_in_chain_gemmi(const gemmi::Chain* chain, bool new_res_no_by_hundreds) {
   std::pair<short int, int> p(0, 1);
   int max_res_no = -9999;

   if (chain) {
      size_t nres = chain->residues.size();
      if (nres > 0) {
         for (size_t ires = nres - 1; ires + 1 > 0; ires--) {
            const gemmi::Residue& residue = chain->residues[static_cast<int>(ires)];
            int seqnum = residue.seqid.num.value;
            if (seqnum > max_res_no) {
               max_res_no = seqnum;
               bool is_het_residue_flag = is_het_residue_gemmi(&residue);
               if (is_het_residue_flag) {
                  p = std::pair<short int, int>(1, residue.seqid.num.value + 1);
               } else {
                  if (new_res_no_by_hundreds) {
                     if (max_res_no < 9999) {
                        int res_no = round_up_by_hundreds(max_res_no + 1);
                        p = std::pair<short int, int>(1, res_no + 1);
                     }
                  } else {
                     if (max_res_no < 9999) {
                        p = std::pair<short int, int>(1, max_res_no + 1);
                     }
                  }
               }
            }
         }
         if (!p.first) {
            // first the first space starting from the front
            int test_resno_start = 1001;
            bool is_clear = false;
            while (!is_clear) {
               is_clear = true;
               for (size_t iser = 0; iser < nres; iser++) {
                  int resno_res = chain->residues[static_cast<int>(iser)].seqid.num.value;
                  if (resno_res >= test_resno_start) {
                     if (resno_res <= (test_resno_start + 10)) {
                        is_clear = false;
                     }
                  }
                  if (!is_clear) break;
               }
               test_resno_start += 100;
            }
            p = std::pair<short int, int>(1, test_resno_start);
         }
      }
   }
   return p;
}

} // namespace cho
} // namespace coot