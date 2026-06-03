#include "function.hh"
#include "geometry/protein-geometry.hh"
#include <gemmi/cifdoc.hpp>
#include <cstdlib>

namespace coot {

void add_chem_mod_angle_gemmi(protein_geometry& geom, const gemmi::cif::Loop& loop) {

  for (std::size_t j = 0; j < loop.length(); j++) {

    int ierr_tot = 0;

    std::string mod_id;
    std::string function;
    std::string atom_id_1;
    std::string atom_id_2;
    std::string atom_id_3;
    double new_value_angle = 0.0;
    double new_value_angle_esd = 0.0;

    auto get_str = [&](const std::string& col_name) -> std::string {
        int idx = loop.find_tag(col_name);
        if (idx >= 0) {
            return loop.val(j, static_cast<std::size_t>(idx));
        }
        return "";
    };

    mod_id     = get_str("mod_id");
    function   = get_str("function");
    atom_id_1  = get_str("atom_id_1");
    atom_id_2  = get_str("atom_id_2");
    atom_id_3  = get_str("atom_id_3");

    // Parse new_value_angle
    {
        std::string val = get_str("new_value_angle");
        if (!val.empty() && val != ".") {
            try {
                new_value_angle = std::stod(val);
            } catch (...) {
                ierr_tot += 1;
            }
        } else {
            ierr_tot += 1;
        }
    }

    // Parse new_value_angle_esd
    {
        std::string val = get_str("new_value_angle_esd");
        if (!val.empty() && val != ".") {
            try {
                new_value_angle_esd = std::stod(val);
            } catch (...) {
                ierr_tot += 1;
            }
        } else {
            ierr_tot += 1;
        }
    }

    if (ierr_tot == 0) {
        coot::chem_mod_angle cma(function,
                                 atom_id_mmdb_expand(atom_id_1),
                                 atom_id_mmdb_expand(atom_id_2),
                                 atom_id_mmdb_expand(atom_id_3),
                                 new_value_angle,
                                 new_value_angle_esd);
        geom.mods[mod_id].add_mod_angle(cma);
    }
  }
}

} // namespace coot