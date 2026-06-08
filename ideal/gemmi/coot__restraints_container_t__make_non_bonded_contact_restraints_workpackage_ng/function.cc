#include "function.hh"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace gemmi;

// Inline is_1_4 logic from raic.is_1_4()
// From simple-restraint.cc: checks if i and j share an angle restraint
// by looking at angle_atom_indices[i] and angle_atom_indices[j]
static bool is_1_4_inline(
    int i, int j,
    const std::vector<bool> &fixed_atom_flags,
    const std::vector<int> &angle_atom_indices) {
   bool is_1_4_related = true;
   bool is_1_4_related_via_non_fixed_atom = false;
   bool angle_exists_in_angle_restraints = false;

   if (angle_atom_indices[i] != -1 && angle_atom_indices[j] != -1) {
      angle_exists_in_angle_restraints = true;
   } else if (!fixed_atom_flags[0] && angle_atom_indices[i] == -1) {
      is_1_4_related = false;
   } else if (!fixed_atom_flags[1] && angle_atom_indices[j] == -1) {
      is_1_4_related = false;
   } else if (!angle_exists_in_angle_restraints) {
      is_1_4_related = false;
   }

   return is_1_4_related;
}

void coot::make_non_bonded_contact_restraints_workpackage_ng_gemmi(
    int ithread,
    int imol,
    const coot::protein_geometry &geom,
    const std::vector<std::set<int> > &bonded_atom_indices,
    const std::vector<int> &angle_atom_indices,
    const std::vector<std::set<unsigned int> > &vcontacts,
    std::pair<unsigned int, unsigned int> atom_index_range_pair,
    const std::set<int> &fixed_atom_indices,
    const std::vector<std::string> &energy_type_for_atom,
    bool extended_atom_mode,
    const std::vector<gemmi::Atom *> &atom,
    const std::vector<gemmi::Residue *> &atom_residues,
    const std::vector<gemmi::Chain *> &atom_chains,
    const std::vector<int> &atom_residue_indices,
    const std::vector<bool> &atom_is_metal,
    const std::vector<bool> &atom_is_hydrogen,
    const std::vector<bool> &H_atom_parent_atom_is_donor_vec,
    const std::vector<bool> &atom_is_acceptor_vec,
    std::vector<std::set<int> > *non_bonded_contacts_atom_indices_p,
    std::vector<simple_restraint> *nbc_restraints_fragment_p,
    std::atomic<unsigned int> &done_count) {

   std::vector<std::set<int> > &non_bonded_contacts_atom_indices = *non_bonded_contacts_atom_indices_p;
   std::map<std::string, std::pair<bool, std::vector<std::list<std::string> > > > residue_ring_map_cache;

   auto tweak_neighbours_for_rn_diff_is_1 = [] (const std::vector<bool> &fixed_atom_flags,
                                                const std::string &atom_name_1,
                                                const std::string &atom_name_2,
                                                const int &res_no_1,
                                                const int &res_no_2,
                                                bool second_is_pro,
                                                double *dist_min_p) {
                                               double &dist_min(*dist_min_p);
                                               bool strange_exception = false;
                                               if (fixed_atom_flags.size()) {
                                                  if (fixed_atom_flags[0] || fixed_atom_flags[1]) {
                                                     if (atom_name_1 == " O  ")
                                                        if (atom_name_2 == " CA ")
                                                           strange_exception = true;
                                                     if (atom_name_1 == " CA ")
                                                        if (atom_name_2 == " O  ")
                                                           strange_exception = true;
                                                     if (atom_name_1 == " N  ")
                                                        if (atom_name_2 == " CB ")
                                                           strange_exception = true;
                                                     if (atom_name_1 == " CB ")
                                                        if (atom_name_2 == " N  ")
                                                           strange_exception = true;
                                                     if (atom_name_1 == " C  ")
                                                        if (atom_name_2 == " CB ")
                                                           strange_exception = true;
                                                  }
                                               }
                                               if (strange_exception)
                                                  dist_min = 2.7;

                                               if (res_no_2 > res_no_1) {
                                                  if (atom_name_1 == " C  ") {
                                                     if (atom_name_2 == " HA " || atom_name_2 == "HA2" || atom_name_2 == " HA3") {
                                                        strange_exception = true;
                                                        dist_min = 2.4;
                                                     }
                                                  }
                                                  if (atom_name_1 == " HA " || atom_name_1 == "HA2" || atom_name_1 == " HA3") {
                                                     if (atom_name_2 == " N  ") {
                                                        strange_exception = true;
                                                        dist_min = 2.41;
                                                     }
                                                  }

                                                  if (atom_name_2 == " C  ") {
                                                     if (atom_name_1 == " H1 " || atom_name_1 == " H2 " || atom_name_1 == " H3 ") {
                                                        dist_min = 2.15;
                                                     }
                                                  }

                                                  if (atom_name_1 == " N  ") {
                                                     if (atom_name_2 == " H  ") {
                                                        strange_exception = true;
                                                        dist_min = 2.4;
                                                     }
                                                  }

                                                  if (second_is_pro) {
                                                     if (atom_name_1 == " CA ") {
                                                        if (atom_name_2 == " CD ") {
                                                           dist_min = 3.1;
                                                        }
                                                     }
                                                  }
                                               } else {
                                                  if (atom_name_1 == " HA " || atom_name_1 == "HA2" || atom_name_1 == " HA3") {
                                                     if (atom_name_2 == " C  ") {
                                                        strange_exception = true;
                                                        dist_min = 2.4;
                                                     }
                                                  }
                                                  if (atom_name_1 == " N  ") {
                                                     if (atom_name_2 == " HA " || atom_name_2 == "HA2" || atom_name_2 == " HA3") {
                                                        strange_exception = true;
                                                        dist_min = 2.41;
                                                     }
                                                  }
                                                  if (atom_name_2 == " N  ") {
                                                     if (atom_name_1 == " H  ") {
                                                        strange_exception = true;
                                                        dist_min = 2.4;
                                                     }
                                                  }
                                               }
                                            };

   auto tweak_neighbours_for_rn_diff_is_2 = [] (const std::vector<bool> &fixed_atom_flags,
                                                const std::string &atom_name_1,
                                                const std::string &atom_name_2,
                                                double *dist_min_p) {
                                               double &dist_min(*dist_min_p);
                                               bool strange_exception = false;
                                               if (fixed_atom_flags.size()) {
                                                  if (fixed_atom_flags[0] || fixed_atom_flags[1]) {
                                                     if (atom_name_1 == " C  ")
                                                        if (atom_name_2 == " N  ")
                                                           strange_exception = true;
                                                     if (atom_name_1 == " N  ")
                                                        if (atom_name_2 == " C  ")
                                                           strange_exception = true;

                                                     if (strange_exception)
                                                        dist_min = 2.7;
                                                  }
                                               }
                                               return dist_min;
                                            };

   for (unsigned int i=atom_index_range_pair.first; i<atom_index_range_pair.second; i++) {
      gemmi::Atom *at_1 = atom[i];
      if (!at_1) {
         std::cout << "ERROR:: make_non_bonded_contact_restraints_workpackage_ng()"
                   << " null atom at index " << i << " in range " << atom_index_range_pair.first
                   << " " << atom_index_range_pair.second << " with n_atoms (bonded_atom_indices size()) "
                   << bonded_atom_indices.size() << std::endl;
         continue;
      }

      // gemmi: skip TER-like atoms (empty name)
      if (at_1->name.empty()) continue;

      const std::set<unsigned int> &n_set = vcontacts[i];
      std::string alt_conf_1;
      if (at_1->altloc != '\0')
         alt_conf_1 = std::string(1, at_1->altloc);

      for (unsigned int j_u : n_set) {
         int j = static_cast<int>(j_u);
         if (j < static_cast<int>(i))
            continue;

         if (bonded_atom_indices[i].find(j) != bonded_atom_indices[i].end())
            continue;

         if (non_bonded_contacts_atom_indices[i].find(j) != non_bonded_contacts_atom_indices[i].end())
            continue;

         gemmi::Atom *at_2 = atom[j];
         if (!at_2) continue;

         std::string alt_conf_2;
         if (at_2->altloc != '\0')
            alt_conf_2 = std::string(1, at_2->altloc);

         // alt conf check
         if (!alt_conf_1.empty())
            if (!alt_conf_2.empty())
               if (alt_conf_1 != alt_conf_2)
                  continue;

         // skip if both fixed
         if (fixed_atom_indices.find(static_cast<int>(i)) != fixed_atom_indices.end())
            if (fixed_atom_indices.find(j) != fixed_atom_indices.end())
               continue;

         gemmi::Residue* res_1 = atom_residues[i];
         gemmi::Residue* res_2 = atom_residues[j];

         std::string res_name_1 = res_1->name;
         std::string res_name_2 = res_2->name;
         int res_no_1 = static_cast<int>(res_1->seqid.num.value);
         int res_no_2 = static_cast<int>(res_2->seqid.num.value);

         bool second_is_pro = false;
         if (res_name_2 == "PRO") second_is_pro = true;
         if (res_name_2 == "HYP") second_is_pro = true;

         std::string element_1 = at_1->element.name();
         std::string element_2 = at_2->element.name();
         const std::string &type_1 = energy_type_for_atom[i];
         const std::string &type_2 = energy_type_for_atom[j];

         std::vector<bool> fixed_atom_flags(2, false);
         if (fixed_atom_indices.find(static_cast<int>(i)) != fixed_atom_indices.end()) fixed_atom_flags[0] = true;
         if (fixed_atom_indices.find(j) != fixed_atom_indices.end()) fixed_atom_flags[1] = true;

         double dist_min = 3.4;

         bool in_same_residue_flag = (res_1 == res_2);
         bool in_same_ring_flag = true;

         if (res_2 != res_1) {
            in_same_ring_flag = false;
            in_same_residue_flag = false;
         }

         std::string atom_name_1 = at_1->name;
         std::string atom_name_2 = at_2->name;

         if (in_same_ring_flag)
            in_same_ring_flag = coot::is_in_same_ring_gemmi(imol, *res_2, residue_ring_map_cache,
                                                            atom_name_1, atom_name_2, geom);

         bool is_1_4_related = is_1_4_inline(i, j, fixed_atom_flags, angle_atom_indices);

         bool mc_atoms_tandem = false;
         bool mc_CC_atoms_tandem = false;
         if (!is_1_4_related) {
            int ri_1 = atom_residue_indices[i];
            int ri_2 = atom_residue_indices[j];

            if (atom_name_1 == " C  ")
               if (atom_name_2 == " C  ")
                  if (ri_2 - ri_1 == 1) {
                     mc_atoms_tandem = true;
                     mc_CC_atoms_tandem = true;
                  }
            if (atom_name_1 == " N  ")
               if (atom_name_2 == " N  ")
                  if (ri_2 - ri_1 == -1)
                     mc_atoms_tandem = true;
            if (atom_name_1 == " CA ")
               if (atom_name_2 == " CA ")
                  mc_atoms_tandem = true;
         }

         if (mc_atoms_tandem)
            is_1_4_related = true;

         if (is_1_4_related) {
            if (mc_atoms_tandem) dist_min = 2.99;

            if (in_same_ring_flag) {
               dist_min = 2.64;
            } else {
               if (mc_CC_atoms_tandem) {
                  dist_min = 3.05;
               } else {
                  dist_min = 2.82;
               }
            }

            if (atom_is_hydrogen[i] && atom_is_hydrogen[j]) {
               dist_min = 2.42;
            } else {
               if (atom_is_hydrogen[i])
                  dist_min -= 0.6;
               if (atom_is_hydrogen[j])
                  dist_min -= 0.6;
            }

         } else {
            std::pair<bool, double> nbc_dist = geom.get_nbc_dist_v2(type_1, type_2,
                                                                    element_1, element_2,
                                                                    atom_is_metal[i],
                                                                    atom_is_metal[j],
                                                                    extended_atom_mode,
                                                                    in_same_residue_flag,
                                                                    in_same_ring_flag);

            if (nbc_dist.first) {
               gemmi::CRA cra_1{atom_chains[i], atom_residues[i], at_1};
               gemmi::CRA cra_2{atom_chains[j], atom_residues[j], at_2};

               bool is_O_C_1_5_related = coot::check_for_O_C_1_5_relation_gemmi(cra_1, cra_2);

               if (is_O_C_1_5_related) {
                  dist_min = 2.84;
               } else {
                  dist_min = nbc_dist.second;

                  int rn_diff = std::abs(res_no_2 - res_no_1);

                  if (rn_diff == 1) {
                     tweak_neighbours_for_rn_diff_is_1(fixed_atom_flags, atom_name_1, atom_name_2, res_no_1, res_no_2,
                                                       second_is_pro, &dist_min);
                  }
                  if (rn_diff == 2) {
                     tweak_neighbours_for_rn_diff_is_2(fixed_atom_flags, atom_name_1, atom_name_2, &dist_min);
                  }
               }
            } else {
               dist_min = 2.8;
            }
         }

         bool is_H_non_bonded_contact = false;

         if (atom_is_hydrogen[i]) {
            is_H_non_bonded_contact = true;
            if (H_atom_parent_atom_is_donor_vec[i])
               if (atom_is_acceptor_vec[j])
                  dist_min -= 0.7;
         }

         if (atom_is_hydrogen[j]) {
            is_H_non_bonded_contact = true;
            if (H_atom_parent_atom_is_donor_vec[j])
               if (atom_is_acceptor_vec[i])
                  dist_min -= 0.7;
         }

         if (atom_is_hydrogen[i])
            if (atom_is_hydrogen[j])
               dist_min -= 0.2;

         non_bonded_contacts_atom_indices[i].insert(j);
         simple_restraint::nbc_function_t nbcf = simple_restraint::LENNARD_JONES;
         simple_restraint r(NON_BONDED_CONTACT_RESTRAINT,
                            nbcf, i, j,
                            is_H_non_bonded_contact,
                            fixed_atom_flags, dist_min);
         nbc_restraints_fragment_p->push_back(r);
      }
   }

   done_count += 1;
}
