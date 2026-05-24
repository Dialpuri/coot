#include "function.hh"

namespace coot {

static bool is_empty_ins_code(const std::string& ins_code) {
  return ins_code.empty() || ins_code == " ";
}

std::vector<cablam_markup_t>
make_cablam_markups_gemmi(const std::vector<std::pair<residue_spec_t, double> >& residues,
                          gemmi::Structure& st) {
  std::vector<cablam_markup_t> v;
  
  for (const auto& pair : residues) {
    const residue_spec_t& spec = pair.first;
    double score = pair.second;
    
    gemmi::Chain* found_chain = nullptr;
    gemmi::Residue* found_res = nullptr;
    
    for (auto& model : st.models) {
      for (auto& chain : model.chains) {
        for (auto& res : chain.residues) {
          std::string res_ins_code = (res.seqid.icode == ' ') ? "" : std::string(1, res.seqid.icode);
          std::string spec_ins_code = is_empty_ins_code(spec.ins_code) ? "" : spec.ins_code;
          
          if (chain.name == spec.chain_id && 
              res.seqid.num.value == spec.res_no && 
              res_ins_code == spec_ins_code) {
            found_chain = &chain;
            found_res = &res;
            break;
          }
        }
        if (found_res) break;
      }
      if (found_res) break;
    }
    
    if (!found_res || !found_chain) {
      continue;
    }
    
    int res_idx = -1;
    for (size_t i = 0; i < found_chain->residues.size(); ++i) {
      if (&found_chain->residues[i] == found_res) {
        res_idx = static_cast<int>(i);
        break;
      }
    }
    
    if (res_idx < 0) continue;
    if (found_chain->residues.size() < 4) continue;
    
    gemmi::CRA cra{found_chain, found_res, nullptr};
    cablam_markup_t cm;
    cm.cra = cra;
    cm.score = score;
    v.push_back(cm);
  }
  
  return v;
}

} // namespace coot