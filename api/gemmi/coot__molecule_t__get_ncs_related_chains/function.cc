#include "function.hh"
#include <map>
#include <string>

namespace coot {
namespace molecule_t {

std::vector<std::vector<const gemmi::Chain*>>
ncs_related_chains_gemmi(const gemmi::Model& model) {
   std::vector<std::vector<const gemmi::Chain*>> v;

   auto chains_match = [] (const gemmi::Chain* c1, const gemmi::Chain* c2) {
      size_t n_res_1 = c1->residues.size();
      size_t n_res_2 = c2->residues.size();
      std::map<int, std::string> m1;
      std::map<int, std::string> m2;
      for (size_t ires=0; ires<n_res_1; ires++) {
         const gemmi::Residue& res = c1->residues[ires];
         std::string res_name = res.name;
         int res_no = res.seqid.num.value;
         if (res_name != "HOH")
            m1[res_no] = res_name;
      }
      for (size_t ires=0; ires<n_res_2; ires++) {
         const gemmi::Residue& res = c2->residues[ires];
         std::string res_name = res.name;
         int res_no = res.seqid.num.value;
         if (res_name != "HOH")
            m2[res_no] = res_name;
      }

      std::map<int, std::string>::const_iterator it_1;
      size_t n_count =  0;
      size_t n_match =  0;
      for (it_1=m1.begin(); it_1!=m1.end(); ++it_1) {
         n_count++;
         std::map<int, std::string>::const_iterator it_2 = m2.find(it_1->first);
         if (it_2 != m2.end()) {
            if (it_2->second == it_1->second)
               n_match++;
         }
      }
      if (n_count > 0) {
         float f1 = static_cast<float>(n_count);
         float f2 = static_cast<float>(n_match);
         if (f2/f1 > 0.7f)
            return true;
      }
      return false;
   };

   for (const auto& chain : model.chains) {
      if (chain.residues.size() > 0) {
         bool found_a_match = false;
         for (size_t iv1=0; iv1<v.size(); ++iv1) {
            std::vector<const gemmi::Chain*>& vv = v[iv1];
            const gemmi::Chain* vv_chain_p = vv[0];
            bool cm = chains_match(vv_chain_p, &chain);
            if (cm) {
               vv.push_back(&chain);
               found_a_match = true;
               break;
            }
         }
         if (! found_a_match) {
            v.push_back({&chain});
         }
      }
   }

   return v;
}

std::vector<std::vector<std::string>>
get_ncs_related_chains_gemmi(const gemmi::Structure& st) {
   std::vector<std::vector<std::string>> v;
   
   if (st.models.empty())
      return v;
   
   const gemmi::Model& model = st.models[0];
   std::vector<std::vector<const gemmi::Chain*>> ncs_related_chains =
       ncs_related_chains_gemmi(model);
   
   for (const auto& vv : ncs_related_chains) {
      std::vector<std::string> vc;
      for (const auto& c : vv) {
         vc.push_back(c->name);
      }
      v.push_back(vc);
   }
   
   return v;
}

}
}