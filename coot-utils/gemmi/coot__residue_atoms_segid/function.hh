#pragma once
#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <stdexcept>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__is_member_p/gemmi/function.hh"

namespace coot {

std::string residue_atoms_segid_gemmi(const gemmi::Residue& residue) {
  std::vector<std::string> seg_ids;
  for (const auto& atom : residue.atoms) {
    std::string seg_id = residue.segment;
    if (seg_ids.empty()) {
      seg_ids.push_back(seg_id);
    } else {
      if (!coot::is_member_p_gemmi(seg_ids, seg_id)) {
        std::string mess = "No consistent segids for residue ";
        mess += std::to_string(residue.seqid.num.value);
        throw std::runtime_error(mess);
      }
    }
  }

  if (seg_ids.empty()) {
    std::string mess = "No segids for residue ";
    mess += std::to_string(residue.seqid.num.value);
    throw std::runtime_error(mess);
  }

  return seg_ids[0];
}

} // namespace coot