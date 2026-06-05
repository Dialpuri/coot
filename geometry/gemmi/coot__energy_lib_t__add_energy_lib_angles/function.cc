#include "function.hh"

void coot::energy_lib_t::add_energy_lib_angles_gemmi(const gemmi::cif::Loop& loop) {

   // CIF tags include the category prefix: _lib_angle.atom_type_1, etc.
   int col1 = loop.find_tag("_lib_angle.atom_type_1");
   int col2 = loop.find_tag("_lib_angle.atom_type_2");
   int col3 = loop.find_tag("_lib_angle.atom_type_3");
   int col_val = loop.find_tag("_lib_angle.value");
   int col_const = loop.find_tag("_lib_angle.const");

   for (int j = 0; j < int(loop.length()); j++) {

      std::string atom_type_1;
      std::string atom_type_2;
      std::string atom_type_3;
      double value = 90.0;
      double value_esd = 1.8;
      double ktheta = 0;
      int ierr;
      int ierr_tot = 0;

      if (col1 >= 0) {
          std::string val = loop.val(j, col1);
          if (!val.empty() && val != ".") atom_type_1 = val;
      }

      if (col2 >= 0) {
          std::string val = loop.val(j, col2);
          if (!val.empty() && val != ".") atom_type_2 = val;
      }

      if (col3 >= 0) {
          std::string val = loop.val(j, col3);
          if (!val.empty() && val != ".") atom_type_3 = val;
      }

      // Get "value" column
      if (col_val >= 0) {
          std::string val_str = loop.val(j, col_val);
          if (!val_str.empty() && val_str != "." && val_str != "?") {
              try {
                  value = std::stod(val_str);
              } catch (...) {
                  // ignore
              }
          }
      }

      // Get "const" column — MMDB GetReal sets to 0 on "." or "?"
      ierr = 0;
      if (col_const >= 0) {
          std::string val_str = loop.val(j, col_const);
          if (!val_str.empty() && val_str != "." && val_str != "?") {
              try {
                  ktheta = std::stod(val_str);
              } catch (...) {
                  ktheta = 0;
                  ierr = 1;
              }
          } else {
              ktheta = 0;
              ierr = 1;
          }
      } else {
          ktheta = 0;
          ierr = 1;
      }

      if (ierr != 0)
         value_esd = 2.8;
      else
         value_esd = ktheta * 0.04;

      if (ierr_tot == 0) {

         energy_lib_angle angle(atom_type_1, atom_type_2, atom_type_3, ktheta, value, value_esd);
         angles.push_back(angle);
      } else {
         std::cout << "  reject energy lib angle "
                   << "\"" << atom_type_1 << "\" "
                   << "\"" << atom_type_2 << "\" "
                   << "\"" << atom_type_3 << "\" "
                   << std::endl;
      }
   }
}