#include "function.hh"

namespace coot::protein_geometry {

void link_torsion_gemmi(const gemmi::cif::Loop &loop) {

   std::string link_id;
   std::string atom_id_1, atom_id_2, atom_id_3, atom_id_4;
   double value_angle, value_angle_esd;
   int atom_1_comp_id, atom_2_comp_id, atom_3_comp_id, atom_4_comp_id;
   int period;
   std::string id("unknown");

   for (size_t j = 0; j < loop.length(); j++) {

      // Get string fields
      int col_link_id = loop.find_tag("link_id");
      if (col_link_id >= 0) link_id = loop.val(j, col_link_id);

      id = "unknown";
      int col_id = loop.find_tag("id");
      if (col_id >= 0) id = loop.val(j, col_id);

      int col_atom_id_1 = loop.find_tag("atom_id_1");
      if (col_atom_id_1 >= 0) atom_id_1 = loop.val(j, col_atom_id_1);

      int col_atom_id_2 = loop.find_tag("atom_id_2");
      if (col_atom_id_2 >= 0) atom_id_2 = loop.val(j, col_atom_id_2);

      int col_atom_id_3 = loop.find_tag("atom_id_3");
      if (col_atom_id_3 >= 0) atom_id_3 = loop.val(j, col_atom_id_3);

      int col_atom_id_4 = loop.find_tag("atom_id_4");
      if (col_atom_id_4 >= 0) atom_id_4 = loop.val(j, col_atom_id_4);

      // Get integer fields
      int col_a1 = loop.find_tag("atom_1_comp_id");
      int col_a2 = loop.find_tag("atom_2_comp_id");
      int col_a3 = loop.find_tag("atom_3_comp_id");
      int col_a4 = loop.find_tag("atom_4_comp_id");
      int col_val = loop.find_tag("value_angle");
      int col_esd = loop.find_tag("value_angle_esd");
      int col_per = loop.find_tag("period");

      bool ok = true;

      if (col_a1 < 0 || col_a2 < 0 || col_a3 < 0 || col_a4 < 0 ||
          col_val < 0 || col_esd < 0 || col_per < 0) {
         ok = false;
      }

      if (ok) {
         atom_1_comp_id = gemmi::cif::as_int(loop.val(j, col_a1));
         atom_2_comp_id = gemmi::cif::as_int(loop.val(j, col_a2));
         atom_3_comp_id = gemmi::cif::as_int(loop.val(j, col_a3));
         atom_4_comp_id = gemmi::cif::as_int(loop.val(j, col_a4));
         value_angle = std::stod(loop.val(j, col_val));
         value_angle_esd = std::stod(loop.val(j, col_esd));
         period = gemmi::cif::as_int(loop.val(j, col_per));
      }

      if (ok) {
         link_add_torsion_gemmi(
             link_id,
             atom_1_comp_id, atom_2_comp_id, atom_3_comp_id, atom_4_comp_id,
             atom_id_1, atom_id_2, atom_id_3, atom_id_4,
             value_angle, value_angle_esd, period, id);
      } else {
         std::cout << "problem reading link torsion mmCIFLoop" << std::endl;
      }
   }
}

} // namespace coot::protein_geometry