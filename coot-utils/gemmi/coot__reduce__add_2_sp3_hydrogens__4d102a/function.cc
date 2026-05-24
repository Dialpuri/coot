#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/math.hpp>
#include <iostream>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_alt_confs/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_pair_by_bisection/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_hydrogen_atom/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {
namespace reduce {

void add_2_sp3_hydrogens_gemmi(
    const std::string &H_at_name_1,
    const std::string &H_at_name_2,
    const std::string &at_name_1,
    const std::string &at_name_2,
    const std::string &at_name_3,
    double bond_length,
    double angle_between_Hs,
    gemmi::Residue *residue_p,
    bool choose_only_farthest_position) {

   std::vector<std::string> alt_confs = coot::util::get_residue_alt_confs_gemmi(*residue_p);
   for (const std::string& alt_conf : alt_confs) {
      gemmi::Atom* at_1 = nullptr;
      gemmi::Atom* at_2 = nullptr;
      gemmi::Atom* at_3 = nullptr;
      
      for (gemmi::Atom& atom : residue_p->atoms) {
         if (atom.name == at_name_1 && atom.altloc == alt_conf[0]) at_1 = &atom;
         if (atom.name == at_name_2 && atom.altloc == alt_conf[0]) at_2 = &atom;
         if (atom.name == at_name_3 && atom.altloc == alt_conf[0]) at_3 = &atom;
      }
      
      if (at_1 && at_2 && at_3) {
         clipper::Coord_orth pos_1(at_1->pos.x, at_1->pos.y, at_1->pos.z);
         clipper::Coord_orth pos_2(at_2->pos.x, at_2->pos.y, at_2->pos.z);
         clipper::Coord_orth pos_3(at_3->pos.x, at_3->pos.y, at_3->pos.z);
         
         std::pair<clipper::Coord_orth, clipper::Coord_orth> Hs =
            position_pair_by_bisection_gemmi(at_1, at_2, at_3, bond_length,
                                       clipper::Util::d2rad(angle_between_Hs));
         
         float bf = at_2->b_iso;
         if (!choose_only_farthest_position) {
            add_hydrogen_atom_gemmi(H_at_name_1, Hs.first, bf, alt_conf, *residue_p);
            add_hydrogen_atom_gemmi(H_at_name_2, Hs.second, bf, alt_conf, *residue_p);
         } else {
            clipper::Coord_orth at_pos_1 = co_gemmi(at_1);
            clipper::Coord_orth at_pos_3 = co_gemmi(at_3);
            clipper::Coord_orth mp(0.5 * (at_pos_1.x() + at_pos_3.x()),
                                   0.5 * (at_pos_1.y() + at_pos_3.y()),
                                   0.5 * (at_pos_1.z() + at_pos_3.z()));
            double d1 = Hs.first.lengthsq();
            double d2 = Hs.second.lengthsq();
            if (d1 > d2)
               add_hydrogen_atom_gemmi(H_at_name_1, Hs.first, bf, alt_conf, *residue_p);
            else
               add_hydrogen_atom_gemmi(H_at_name_1, Hs.second, bf, alt_conf, *residue_p);
         }
      } else {
         if (!alt_conf.empty()) {
            std::cout << "WARNING: in add_2_sp3_hydrogens_gemmi() could not find atoms "
                      << at_name_1 << " " << at_name_2 << " " << at_name_3
                      << " for alt-conf \"" << alt_conf << "\"" << std::endl;
         }
      }
   }
}

void add_2_sp3_hydrogens_gemmi(
    const std::string &H_at_name_1,
    const std::string &H_at_name_2,
    const std::string &first_neighb,
    const std::vector<std::string> &second_neighb_vec,
    double bond_length,
    double angle_between_Hs,
    gemmi::Residue *residue_p,
    bool choose_only_farthest_position) {

   if (second_neighb_vec.size() == 2) {
      const std::string &second_1 = second_neighb_vec[0];
      const std::string &second_2 = second_neighb_vec[1];
      add_2_sp3_hydrogens_gemmi(H_at_name_1, H_at_name_2, second_1, first_neighb, second_2, bond_length,
                          angle_between_Hs, residue_p, choose_only_farthest_position);
   } else {
      std::cout << "WARNING:: in add_2_sp3_hydrogens_gemmi() second_neighb_vec.size() is "
                << second_neighb_vec.size() << std::endl;
   }
}

}
}