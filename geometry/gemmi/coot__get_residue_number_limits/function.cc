#include "function.hh"

namespace coot {

std::map<std::string, std::pair<int, int> >
get_residue_number_limits_gemmi(const gemmi::Structure& st) {

   std::map<std::string, std::pair<int, int> > limits;
   if (!st.models.empty()) {
      const gemmi::Model& model = st.models[0];
      int n_chains = (int)model.chains.size();
      for (int ichain = 0; ichain < n_chains; ichain++) {
         const gemmi::Chain& chain = model.chains[ichain];
         std::string chain_id = chain.name;
         int nres = (int)chain.residues.size();
         int min_res_no =  9999999;
         int max_res_no = -9999999;
         for (int ires = 0; ires < nres; ires++) {
            const gemmi::Residue& residue = chain.residues[ires];
            int res_no = residue.seqid.num.value;
            if (res_no < min_res_no) min_res_no = res_no;
            if (res_no > max_res_no) max_res_no = res_no;
         }
         if (min_res_no < 9999999) {
            if (max_res_no > -9999999) {
               std::pair<int, int> res_lims(min_res_no, max_res_no);
               limits[chain_id] = res_lims;
            }
         }
      }
   }
   return limits;
}

} // namespace coot