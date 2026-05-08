#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {
  struct molecule_t {
    int get_number_of_hydrogen_atoms_gemmi(const gemmi::Structure& st) const {
      if (st.models.empty())
        return -1;
      int n = 0;
      for (const gemmi::Chain& chain : st.models[0].chains)
        for (const gemmi::Residue& res : chain.residues)
          for (const gemmi::Atom& atom : res.atoms)
            if (atom.element.is_hydrogen() && res.name != "TER")
              n++;
      return n;
    }
  };
}