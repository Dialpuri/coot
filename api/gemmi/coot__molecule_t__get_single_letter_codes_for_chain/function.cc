#include "function.hh"
#include "coot/coot-utils/coot-coord-utils.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__sort_residues_by_seqno/gemmi/function.hh"

namespace coot {

std::vector<std::pair<residue_spec_t, std::string>>
get_single_letter_codes_for_chain_gemmi(
  const gemmi::Model &model, const std::string &chain_id) {

  std::vector<std::pair<residue_spec_t, std::string>> v;

  for (const gemmi::Chain &chain : model.chains) {
    if (chain.name == chain_id) {
      std::vector<gemmi::Residue*> residue_ptrs;
      for (const gemmi::Residue &res : chain.residues) {
        residue_ptrs.push_back(const_cast<gemmi::Residue*>(&res));
      }

      auto sorted = coot::util::sort_residues_by_seqno(
        residue_ptrs.data(), static_cast<int>(residue_ptrs.size()));

      for (const auto &pair : sorted) {
        const gemmi::Residue *res = pair.first;
        residue_spec_t res_spec;
        res_spec.chain_id = chain.name;
        res_spec.res_no   = res->seqid.num.value;
        res_spec.ins_code = std::string(1, res->seqid.icode);
        if (res_spec.ins_code == " ") res_spec.ins_code = "";

        std::string res_name = res->name;
        std::string s = coot::util::three_letter_to_one_letter_with_specials(res_name);
        v.push_back({res_spec, s});
      }
    }
  }

  return v;
}

} // namespace coot