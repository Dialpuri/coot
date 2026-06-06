#pragma once

#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"

namespace coot {

class beam_in_linked_residue {
public:
  gemmi::Residue* residue_ref;
  gemmi::Residue* template_res_ref;
  gemmi::Residue* template_res_mov;
  std::string comp_id_ref;
  std::string comp_id_new;
  std::string link_type;
  bool have_template;
  coot::protein_geometry* geom_p;

  beam_in_linked_residue(gemmi::Residue* res, const std::string& link,
                         const std::string& comp_id,
                         coot::protein_geometry* geom);

  bool setup_by_comp_id(const std::string& template_id,
                        const std::string& reference_id);

  gemmi::Residue* get_residue_raw_gemmi() const;
};

} // namespace coot