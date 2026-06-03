#include "function.hh"

void add_chem_mod_tor_gemmi(
    coot::protein_geometry& geom,
    const gemmi::cif::Loop& cif_loop) {

  for (size_t j = 0; j < cif_loop.length(); j++) {

    int ierr_tot = 0;
    int ierr = 0;

    std::string mod_id;
    std::string function;
    std::string atom_id_1;
    std::string atom_id_2;
    std::string atom_id_3;
    std::string atom_id_4;
    double new_value_angle = 0.0;
    double new_value_angle_esd = 0.0;
    int new_period = 0;

    // Try to find column indices
    int mod_id_col = cif_loop.find_tag("mod_id");
    int function_col = cif_loop.find_tag("function");
    int atom_1_col = cif_loop.find_tag("atom_id_1");
    int atom_2_col = cif_loop.find_tag("atom_id_2");
    int atom_3_col = cif_loop.find_tag("atom_id_3");
    int atom_4_col = cif_loop.find_tag("atom_id_4");
    int angle_col = cif_loop.find_tag("new_value_angle");
    int esd_col = cif_loop.find_tag("new_value_angle_esd");
    int period_col = cif_loop.find_tag("new_period");

    // Read strings
    if (mod_id_col >= 0) {
      const std::string& v = cif_loop.val(j, mod_id_col);
      mod_id = gemmi::cif::as_string(v);
      if (gemmi::cif::is_null(v)) ierr++;
    } else {
      ierr++;
    }

    if (function_col >= 0) {
      const std::string& v = cif_loop.val(j, function_col);
      function = gemmi::cif::as_string(v);
      if (gemmi::cif::is_null(v)) ierr++;
    } else {
      ierr++;
    }

    if (atom_1_col >= 0) {
      const std::string& v = cif_loop.val(j, atom_1_col);
      atom_id_1 = gemmi::cif::as_string(v);
      if (gemmi::cif::is_null(v)) ierr++;
    } else {
      ierr++;
    }

    if (atom_2_col >= 0) {
      const std::string& v = cif_loop.val(j, atom_2_col);
      atom_id_2 = gemmi::cif::as_string(v);
      if (gemmi::cif::is_null(v)) ierr++;
    } else {
      ierr++;
    }

    if (atom_3_col >= 0) {
      const std::string& v = cif_loop.val(j, atom_3_col);
      atom_id_3 = gemmi::cif::as_string(v);
      if (gemmi::cif::is_null(v)) ierr++;
    } else {
      ierr++;
    }

    if (atom_4_col >= 0) {
      const std::string& v = cif_loop.val(j, atom_4_col);
      atom_id_4 = gemmi::cif::as_string(v);
      if (gemmi::cif::is_null(v)) ierr++;
    } else {
      ierr++;
    }

    // Read real: new_value_angle
    if (angle_col >= 0) {
      const std::string& v = cif_loop.val(j, angle_col);
      if (gemmi::cif::is_null(v)) {
        ierr++;
      } else {
        new_value_angle = std::stod(v);
      }
    } else {
      ierr++;
    }

    // Read real: new_value_angle_esd
    if (esd_col >= 0) {
      const std::string& v = cif_loop.val(j, esd_col);
      if (gemmi::cif::is_null(v)) {
        ierr++;
      } else {
        new_value_angle_esd = std::stod(v);
      }
    } else {
      ierr++;
    }

    // Read integer: new_period
    if (period_col >= 0) {
      const std::string& v = cif_loop.val(j, period_col);
      if (gemmi::cif::is_null(v)) {
        ierr++;
      } else {
        new_period = gemmi::cif::as_int(v);
      }
    } else {
      ierr++;
    }

    ierr_tot += ierr;

    if (ierr_tot == 0) {
      coot::chem_mod_tor cmt(function,
              coot::atom_id_mmdb_expand(atom_id_1),
              coot::atom_id_mmdb_expand(atom_id_2),
              coot::atom_id_mmdb_expand(atom_id_3),
              coot::atom_id_mmdb_expand(atom_id_4),
              new_value_angle,
              new_value_angle_esd,
              new_period);
      geom.mods[mod_id].add_mod_tor(cmt);
    }
  }
}