#pragma once

#include <gemmi/model.hpp>
#include <string>

namespace coot {

class protein_geometry;
class chem_link;
class dict_link_info_t {
public:
  bool check_for_order_switch_gemmi(const gemmi::Residue *residue_ref,
                                    const gemmi::Residue *residue_new,
                                    const std::string &link_type,
                                    const coot::protein_geometry &geom) const;
};

} // namespace coot