#include "function.hh"

#include "coot/utils/coot-utils.hh"

#include <cmath>

namespace coot {

std::string format_gemmi(
    simple_restraint const& r,
    std::vector<atom_spec_t> const& atoms,
    double distortion) {

   std::string s;

   if (r.restraint_type == BOND_RESTRAINT) {
      s = "Bond ";
      atom_spec_t spec_1 = atoms[r.atom_index_1];
      atom_spec_t spec_2 = atoms[r.atom_index_2];
      s += spec_1.label();
      s += " ";
      s += spec_2.label();
      s += " ";
      if (distortion >= 0) {
         s += "  z = ";
         s += util::float_to_string_using_dec_pl(std::sqrt(distortion), 2);
      }
   }
   if (r.restraint_type == ANGLE_RESTRAINT) {
      s = "Angle ";
      atom_spec_t spec_1 = atoms[r.atom_index_1];
      atom_spec_t spec_2 = atoms[r.atom_index_2];
      atom_spec_t spec_3 = atoms[r.atom_index_3];
      s += spec_1.label();
      s += " ";
      s += spec_2.label();
      s += " ";
      s += spec_3.label();
      s += " ";
      if (distortion >= 0) {
         s += "  z = ";
         s += util::remove_whitespace(util::float_to_string_using_dec_pl(std::sqrt(distortion), 2));
      }
   }
   if (r.restraint_type == TORSION_RESTRAINT) {
      s = "Torsion ";
      atom_spec_t spec_1 = atoms[r.atom_index_1];
      atom_spec_t spec_2 = atoms[r.atom_index_2];
      atom_spec_t spec_3 = atoms[r.atom_index_3];
      atom_spec_t spec_4 = atoms[r.atom_index_4];
      s += spec_1.label();
      s += " ";
      s += spec_2.label();
      s += " ";
      s += spec_3.label();
      s += " ";
      s += spec_4.label();
      s += " ";
      if (distortion >= 0) {
         s += "  z = ";
         s += util::remove_whitespace(util::float_to_string_using_dec_pl(std::sqrt(distortion), 2));
      }
   }
   if (r.restraint_type == TRANS_PEPTIDE_RESTRAINT) {
      s = "Trans_peptide ";
      atom_spec_t spec_1 = atoms[r.atom_index_1];
      atom_spec_t spec_2 = atoms[r.atom_index_2];
      atom_spec_t spec_3 = atoms[r.atom_index_3];
      atom_spec_t spec_4 = atoms[r.atom_index_4];
      s += spec_1.label();
      s += " ";
      s += spec_2.label();
      s += " ";
      s += spec_3.label();
      s += " ";
      s += spec_4.label();
      s += " ";
      if (distortion >= 0) {
         s += "  score: ";
         s += util::remove_whitespace(util::float_to_string_using_dec_pl(distortion, 2));
         s += " (non-sqrt)";
      }
   }
   if (r.restraint_type == PLANE_RESTRAINT) {
      s = "Plane ";
      for (std::size_t j = 0; j < r.plane_atom_index.size(); j++) {
         atom_spec_t spec = atoms[r.plane_atom_index[j].first];
         s += spec.label();
         s += " ";
      }
      if (distortion >= 0) {
         s += "  z = ";
         s += util::remove_whitespace(util::float_to_string_using_dec_pl(std::sqrt(distortion), 2));
      }
   }
   if (r.restraint_type == NON_BONDED_CONTACT_RESTRAINT) {
      s = "Non-Bonded-Contact ";
      atom_spec_t spec_1 = atoms[r.atom_index_1];
      atom_spec_t spec_2 = atoms[r.atom_index_2];
      s += spec_1.label();
      s += " ";
      s += spec_2.label();
      s += " ";
      if (distortion >= 0) {
         s += "  z = ";
         s += util::remove_whitespace(util::float_to_string_using_dec_pl(std::sqrt(distortion), 2));
      }
   }
   if (r.restraint_type == CHIRAL_VOLUME_RESTRAINT) {
      s = "Chiral ";
      atom_spec_t spec_c = atoms[r.atom_index_centre];
      atom_spec_t spec_1 = atoms[r.atom_index_1];
      atom_spec_t spec_2 = atoms[r.atom_index_2];
      atom_spec_t spec_3 = atoms[r.atom_index_3];
      s += spec_c.label();
      s += " neigbhs: ";
      s += spec_1.label();
      s += " ";
      s += spec_2.label();
      s += " ";
      s += spec_3.label();
      if (distortion >= 0) {
         s += "  z = ";
         s += util::remove_whitespace(util::float_to_string_using_dec_pl(std::sqrt(distortion), 2));
      }
   }
   if (r.restraint_type == RAMACHANDRAN_RESTRAINT) {
      s = "Rama ";
      atom_spec_t spec_1 = atoms[r.atom_index_1];
      atom_spec_t spec_2 = atoms[r.atom_index_2];
      atom_spec_t spec_3 = atoms[r.atom_index_3];
      atom_spec_t spec_4 = atoms[r.atom_index_4];
      atom_spec_t spec_5 = atoms[r.atom_index_5];
      s += spec_1.label();
      s += " ";
      s += spec_2.label();
      s += " ";
      s += spec_3.label();
      s += " ";
      s += spec_4.label();
      s += " ";
      s += spec_5.label();
      s += " ";
      s += util::remove_whitespace(util::float_to_string_using_dec_pl(distortion, 2));
   }

   if (r.restraint_type == TARGET_POS_RESTRAINT) {
      s = "Target_pos ";
      atom_spec_t spec_1 = atoms[r.atom_index_1];
      s += spec_1.label();
      s += " ";
      s += util::remove_whitespace(util::float_to_string_using_dec_pl(distortion, 2));
   }

   return s;
}

} // namespace coot