#include "function.hh"
#include "coot-utils/coot-h-bonds.hh"
#include <iostream>
#include <gemmi/model.hpp>

namespace coot {

double distance_gemmi(const gemmi::Atom* at_1, const gemmi::Atom* at_2);
double angle_gemmi(const gemmi::Atom *at_1, const gemmi::Atom *at_2, const gemmi::Atom *at_3);

namespace h_bonds_gemmi {

std::pair<bool, h_bond_gemmi> 
make_h_bond_from_environment_residue_hydrogen_gemmi(const gemmi::Atom* at_1, // acceptor on ligand
                                                    const gemmi::Atom* at_2, // H on residue
                                                    const std::vector<std::pair<const gemmi::Atom*, float>>& nb_1,
                                                    const std::vector<std::pair<const gemmi::Atom*, float>>& nb_2) {
   bool debug = false;
   // In gemmi, atoms don't have parent pointers, so we use empty strings
   std::string resname_1 = "";
   std::string resname_2 = "";
   
   if (debug)
      std::cout << "\nDEBUG:: start make_h_bond_from_environment_residue_hydrogen() with"
                << " at_1: " << at_1->name << " " << resname_1
                << " at_2: " << at_2->name << " " << resname_2
                << " nb_1.size(): " << nb_1.size() << " nb_2.size() " << nb_2.size()
                << std::endl;

   double water_dist_max = 3.25; // pass this

   if (debug) {
      for (size_t i=0; i<nb_1.size(); i++)
         std::cout << "    nb of at_1: " << nb_1[i].first->name << std::endl;
      for (size_t i=0; i<nb_2.size(); i++)
         std::cout << "    nb of at_2: " << nb_2[i].first->name << std::endl;
   }

   bool ligand_atom_is_H_flag = false;
   h_bond_gemmi bond(at_2, at_1, ligand_atom_is_H_flag); // H atom goes first for this constructor
   bond.dist = distance_gemmi(at_1, at_2);

   bool neighbour_distances_and_angles_are_good = true;
   bool good_donor_acceptor_dist = false;

   // Dist D-A
   //
   for (size_t iD=0; iD<nb_2.size(); iD++) {
      double dist = distance_gemmi(nb_2[iD].first, at_1);
      if (dist < 3.9) { // McDonald and Thornton
         good_donor_acceptor_dist = true;
         break;
      }
   }

   // Dist D-A for HOH acceptor
   //
   // (in this case nb_2.size() is 0)
   //
   if (resname_2 == "HOH") {
      if (bond.dist < water_dist_max) {
         good_donor_acceptor_dist = true;
         bond.donor = at_2;
      }
   }

   // Angle D-H-A
   //
   for (size_t iD=0; iD<nb_2.size(); iD++) { 
      double angle = angle_gemmi(nb_2[iD].first, at_2, at_1);
      if (debug) {
         std::cout << "   H-on-protein angle 1: " << angle << "  ";
         std::cout << " : "
                   << nb_2[iD].first->name << " "
                   << at_2->name << " "
                   << at_1->name << std::endl;
      }
      if (angle < 90) {
         if (debug) std::cout << "DEBUG:: angle-1 bad" << std::endl;
         neighbour_distances_and_angles_are_good = false;
         break;
      } else {
         if (debug) std::cout << "DEBUG:: angle-1 good" << std::endl;
      }
      if (! bond.donor) {
         bond.donor = nb_2[iD].first;
         bond.angle_1 = angle;
      } 
   }

   // Angle H-A-AA
   // 
   bool found_a_goodie_angle_2 = false;
   for (size_t iA=0; iA<nb_1.size(); iA++) { 
      double angle = angle_gemmi(at_2, at_1, nb_1[iA].first);
      if (debug) {
         std::cout << "   H-on-protein angle 2: " << angle << "  ";
         std::cout << " : "
                   << at_2->name << " "
                   << at_1->name << " "
                   << nb_1[iA].first->name << std::endl;
      }
      if (angle < 90) {
         if (debug) std::cout << "DEBUG:: this angle-2 bad" << std::endl;
      } else {
         found_a_goodie_angle_2 = true;
         if (debug) std::cout << "DEBUG:: angle-2 good" << std::endl;
      }
      if (found_a_goodie_angle_2) {
         if (! bond.acceptor) {
            bond.acceptor = at_1;
            bond.angle_2 = angle;
         }
      }
   }
   if (! found_a_goodie_angle_2)
      neighbour_distances_and_angles_are_good = false;

   // Angle D-A-AA
   //
   bool found_a_goodie_angle_3 = false;
   if (nb_2.size() > 0) {
      for (size_t iD=0; iD<nb_2.size(); iD++) {
         for (size_t iA=0; iA<nb_1.size(); iA++) {
            double angle = angle_gemmi(nb_2[iD].first, at_1, nb_1[iA].first);
            if (debug) {
               std::cout << "   H-on-protein angle 3: " << angle << "  ";
               std::cout << " : "
                         << nb_2[iD].first->name << " "
                         << at_1->name << " "
                         << nb_1[iA].first->name << std::endl;
            }
            if (angle < 90) {
               if (debug) std::cout << "DEBUG:: this angle-3 bad" << std::endl;
            } else {
               if (debug) std::cout << "DEBUG:: this angle-3 good" << std::endl;
               found_a_goodie_angle_3 = true;
            }
            if (found_a_goodie_angle_3) {
               if (! bond.acceptor_neigh) {
                  bond.acceptor_neigh = nb_1[iA].first;
                  bond.angle_3 = angle;
               }
            }
         }
      }
      if (! found_a_goodie_angle_3)
         neighbour_distances_and_angles_are_good = false;
   } else {
      // for HOH, there are ne neighbours of the donor atom (nb_2.size() == 0).
      if (nb_1.size() > 0) {
         double angle = angle_gemmi(at_2, at_1, nb_1[0].first);
         if (! bond.acceptor_neigh) {
            bond.acceptor_neigh = nb_1[0].first;
            bond.angle_3 = angle;
         }
      }
   }

   if (debug)
      std::cout << "DEBUG:: in make_h_bond_from_environment_residue_hydrogen() neighbour_distances_and_angles_are_good: "
                << neighbour_distances_and_angles_are_good << " good_donor_acceptor_dist: " << good_donor_acceptor_dist
                << std::endl;

   return std::pair<bool, h_bond_gemmi> (neighbour_distances_and_angles_are_good && good_donor_acceptor_dist, bond);
}

} // namespace h_bonds_gemmi
} // namespace coot