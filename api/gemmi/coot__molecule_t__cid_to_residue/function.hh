#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t_detail {
  // Parse CID of form "//<chain_id>/<resseq>" and find the residue
  // Returns nullptr if not found or if cid is invalid
  inline gemmi::Residue* cid_to_residue_impl(gemmi::Model& model, const std::string& cid) {
    // CID format: "//<chain_id>/<resseq>" (e.g., "//A/10")
    if (cid.size() < 4 || cid[0] != '/' || cid[1] != '/')
      return nullptr;
    
    // Find the position of the second '/' after chain_id
    size_t slash_pos = cid.find('/', 2);
    if (slash_pos == std::string::npos)
      return nullptr;
    
    std::string chain_id = cid.substr(2, slash_pos - 2);
    std::string resseq_str = cid.substr(slash_pos + 1);
    
    // Parse residue sequence number
    int resseq;
    try {
      resseq = std::stoi(resseq_str);
    } catch (...) {
      return nullptr;
    }
    
    // Find the residue in the model
    for (auto& chain : model.chains) {
      if (chain.name == chain_id) {
        for (auto& residue : chain.residues) {
          if (residue.seqid.num.value == resseq) {
            return &residue;
          }
        }
      }
    }
    
    return nullptr;
  }
}

inline gemmi::Residue* cid_to_residue_gemmi(const std::string& cid, gemmi::Model& model) {
  return molecule_t_detail::cid_to_residue_impl(model, cid);
}
}