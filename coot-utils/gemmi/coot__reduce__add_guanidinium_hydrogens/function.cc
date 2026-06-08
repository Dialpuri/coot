#include "function.hh"
#include <clipper/core/clipper_util.h>

namespace coot {
namespace reduce {

// Helper: find atom in gemmi::Residue by name, with optional altloc filtering
// mimics mmdb::Residue::GetAtom(name, 0, altconf)
static gemmi::Atom* find_atom_in_residue(const gemmi::Residue& residue,
                                          const std::string& name,
                                          const std::string& altconf) {
    char altloc = altconf.empty() ? '\0' : altconf[0];
    for (const gemmi::Atom& a : residue.atoms) {
        // gemmi atom name is padded to 4 chars — compare after padding
        std::string a_name = a.name;
        if (a_name == name) {
            if (altconf.empty()) {
                // MMDB: empty alt_conf matches any altloc
                return const_cast<gemmi::Atom*>(&a);
            } else {
                if (a.altloc == altloc)
                    return const_cast<gemmi::Atom*>(&a);
            }
        }
    }
    return nullptr;
}

void add_guanidinium_hydrogens_gemmi(gemmi::Residue& residue) {

   std::vector<std::string> alt_confs = coot::util::get_residue_alt_confs_gemmi(residue);
   for (unsigned int i = 0; i < alt_confs.size(); i++) {

      // HE
      std::string H_at_name = " HE ";
      double bl = 0.86;

      gemmi::Atom* at_n_1 = find_atom_in_residue(residue, " CD ", alt_confs[i]);
      gemmi::Atom* at_n_2 = find_atom_in_residue(residue, " NE ", alt_confs[i]);
      gemmi::Atom* at_n_3 = find_atom_in_residue(residue, " CZ ", alt_confs[i]);
      if (at_n_1 && at_n_2 && at_n_3) {
         double bf = at_n_2->b_iso;
         clipper::Coord_orth H_pos = position_by_bisection_gemmi(at_n_1, at_n_2, at_n_3, bl);
         add_hydrogen_atom_gemmi(H_at_name, H_pos, static_cast<float>(bf), alt_confs[i], residue);
      } else {
         std::cout << "Fail Residue " << residue.name
                   << " alt-conf \"" << alt_confs[i] << "\""
                   << " failed in add_guanidinium_hydrogens " << std::endl;
      }

      // HH[12][12]
      at_n_1 = find_atom_in_residue(residue, " NE ", alt_confs[i]);
      at_n_2 = find_atom_in_residue(residue, " CZ ", alt_confs[i]);
      gemmi::Atom* at_nh1 = find_atom_in_residue(residue, " NH1", alt_confs[i]);
      gemmi::Atom* at_nh2 = find_atom_in_residue(residue, " NH2", alt_confs[i]);
      if (at_n_1 && at_n_2 && at_nh1 && at_nh2) {
         double bf_nh1 = at_nh1->b_iso;
         double bf_nh2 = at_nh2->b_iso;
         double a = clipper::Util::d2rad(120);
         double t = clipper::Util::d2rad(180);
         clipper::Coord_orth hh11 = position_by_bond_length_angle_torsion_gemmi(at_n_1, at_n_2, at_nh1, bl, a, 0);
         clipper::Coord_orth hh12 = position_by_bond_length_angle_torsion_gemmi(at_n_1, at_n_2, at_nh1, bl, a, t);
         clipper::Coord_orth hh21 = position_by_bond_length_angle_torsion_gemmi(at_n_1, at_n_2, at_nh2, bl, a, 0);
         clipper::Coord_orth hh22 = position_by_bond_length_angle_torsion_gemmi(at_n_1, at_n_2, at_nh2, bl, a, t);
         add_hydrogen_atom_gemmi("HH11", hh11, static_cast<float>(bf_nh1), alt_confs[i], residue);
         add_hydrogen_atom_gemmi("HH12", hh12, static_cast<float>(bf_nh2), alt_confs[i], residue);
         add_hydrogen_atom_gemmi("HH21", hh21, static_cast<float>(bf_nh2), alt_confs[i], residue);
         add_hydrogen_atom_gemmi("HH22", hh22, static_cast<float>(bf_nh2), alt_confs[i], residue);
      }
   }
}

} // namespace reduce
} // namespace coot
