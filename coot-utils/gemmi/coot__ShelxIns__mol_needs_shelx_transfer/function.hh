#pragma once

#include <gemmi/model.hpp>
#include <string>

namespace coot {

class ShelxIns {
public:
  bool mol_needs_shelx_transfer_gemmi(const gemmi::Structure& st) const {
    bool needs_unshelx = true;
    
    // gemmi Structure has models as a vector (0-indexed)
    if (st.models.empty()) {
      return true;
    }
    
    const gemmi::Model& model = st.models[0];
    if (model.chains.size() > 1) {
      needs_unshelx = false;
    } else {
      for (const gemmi::Chain& chain : model.chains) {
        if (!chain.name.empty()) {
          needs_unshelx = false;
          break;
        }
      }
    }
    
    return needs_unshelx;
  }
};

}