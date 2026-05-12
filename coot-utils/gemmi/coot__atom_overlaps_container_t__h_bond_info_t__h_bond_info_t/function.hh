#pragma once

#include <string>

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {

// hb_t enum matching coot/geometry/hb-types.hh
enum hb_t { HB_UNASSIGNED=-1, HB_NEITHER=0, HB_DONOR=1, HB_ACCEPTOR=2, HB_BOTH=3, HB_HYDROGEN=4 };

class atom_overlaps_container_t {
public:
   struct h_bond_info_t {
      bool is_h_bond_H_and_acceptor;
      bool is_h_bond_donor_and_acceptor;
      bool H_is_first_atom_flag;
      bool H_is_second_atom_flag;
      bool donor_is_second_atom_flag;

      h_bond_info_t(int hb_1, int hb_2,
                    const std::string& resname_1,
                    const std::string& resname_2);

      std::string format() const;
   };
};

} // namespace coot