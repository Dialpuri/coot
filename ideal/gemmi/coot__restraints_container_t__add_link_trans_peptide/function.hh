#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <vector>
#include <string>
#include <cmath>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__restraints_container_t__get_asc_index/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

// restraint type enum
enum restraint_type_t {
   ANGLE_RESTRAINT = 0,
   BOND_RESTRAINT,
   PLANAR_RESTRAINT,
   TRANS_PEPTIDE_RESTRAINT,
   NON_BOND_BUMP_RESTRAINT,
   NON_BOND_CONTACT_RESTRAINT,
   HYDROGEN_BOND_RESTRAINT,
   CHIRAL_RESTRAINT,
   NON_BOND_VDW_RESTRAINT,
   NON_BOND_1_4_RESTRAINT,
   DISULFIDE_BOND_RESTRAINT,
   RESTRAINT_TYPE_END
};

class simple_restraint {
public:
   restraint_type_t rest_type;
   int atom_1;
   int atom_2;
   int atom_3;
   int atom_4;
   std::vector<bool> fixed_atom_flags;
   float target;
   float sigma;
   float torsion_restraint_weight;
   int periodicity;

   simple_restraint(restraint_type_t rt, int a1, int a2, int a3, int a4,
                    const std::vector<bool>& f, float t, float s, float w, int p)
      : rest_type(rt), atom_1(a1), atom_2(a2), atom_3(a3), atom_4(a4),
        fixed_atom_flags(f), target(t), sigma(s), torsion_restraint_weight(w),
        periodicity(p) {}
};

class restraints_container_t {
public:

   std::vector<simple_restraint> restraints_vec;

   /**
    * Add a trans peptide (omega) torsion restraint between two residues.
    * Looks for CA-C-N-CA atom pattern (first->CA, first->C, second->N, second->CA).
    * Adds restraint only if C-N distance < 2.0 and torsion is approximately trans (|omega| > pi/2).
    * If add_even_if_cis is true, ignores the torsion check.
    */
   int add_link_trans_peptide_gemmi(
      gemmi::Chain* first_chain,
      gemmi::Residue* first,
      gemmi::Chain* second_chain,
      gemmi::Residue* second,
      bool is_fixed_first,
      bool is_fixed_second,
      bool add_even_if_cis,
      const std::vector<gemmi::CRA>& cras);

   void thread_pool(void*, int) {}

private:
   std::vector<bool> make_fixed_flags(int, int, int, int) const;
   bool add(restraint_type_t, int, int, int, int,
            const std::vector<bool>&, float, float, float, int);
};

inline std::vector<bool> restraints_container_t::make_fixed_flags(
   int index1, int index2, int index3, int index4) const {
   (void)index1; (void)index2; (void)index3; (void)index4;
   return {false, false, false, false};
}

inline bool restraints_container_t::add(
   restraint_type_t rest_type, int atom_1, int atom_2,
   int atom_3, int atom_4,
   const std::vector<bool>& fixed_atom_flags,
   float tar, float sig, float torsion_restraint_weight, int periodicity) {
   bool r = false;
   if (sig > 0.0) {
      restraints_vec.push_back(simple_restraint(rest_type,
                                                atom_1, atom_2, atom_3, atom_4,
                                                fixed_atom_flags, tar, sig,
                                                torsion_restraint_weight, periodicity));
      r = true;
   }
   return r;
}

inline int restraints_container_t::add_link_trans_peptide_gemmi(
   gemmi::Chain* first_chain,
   gemmi::Residue* first,
   gemmi::Chain* second_chain,
   gemmi::Residue* second,
   bool is_fixed_first,
   bool is_fixed_second,
   bool add_even_if_cis,
   const std::vector<gemmi::CRA>& cras) {

   int n_trans_peptide_torsion = 0;

   if (static_cast<int>(first->atoms.size()) <= 0) {
      // no atoms in first residue
   }
   if (static_cast<int>(second->atoms.size()) <= 0) {
      // no atoms in second residue
   }

   std::vector<bool> fixed_flags(4);
   fixed_flags[0] = is_fixed_first;
   fixed_flags[1] = is_fixed_first;
   fixed_flags[2] = is_fixed_second;
   fixed_flags[3] = is_fixed_second;

   // Find CA in first residue
   for (const auto& atom_1 : first->atoms) {
      if (atom_1.name == "CA") {
         // Find C in first residue
         for (const auto& atom_2 : first->atoms) {
            if (atom_2.name == "C") {
               // Find N in second residue
               for (const auto& atom_3 : second->atoms) {
                  if (atom_3.name == "N") {
                     // Find CA in second residue
                     for (const auto& atom_4 : second->atoms) {
                        if (atom_4.name == "CA") {

                           // Build atom specs for index lookup
                           coot::atom_spec_t spec1{atom_1.name, "", first->seqid.num.value, "", first_chain->name};
                           coot::atom_spec_t spec2{atom_2.name, "", first->seqid.num.value, "", first_chain->name};
                           coot::atom_spec_t spec3{atom_3.name, "", second->seqid.num.value, "", second_chain->name};
                           coot::atom_spec_t spec4{atom_4.name, "", second->seqid.num.value, "", second_chain->name};

                           int index1 = coot::get_asc_index_gemmi(cras, spec1);
                           int index2 = coot::get_asc_index_gemmi(cras, spec2);
                           int index3 = coot::get_asc_index_gemmi(cras, spec3);
                           int index4 = coot::get_asc_index_gemmi(cras, spec4);

                           // Get coordinates using gemmi port
                           clipper::Coord_orth a1 = coot::co_gemmi(&atom_1);
                           clipper::Coord_orth a2 = coot::co_gemmi(&atom_2);
                           clipper::Coord_orth a3 = coot::co_gemmi(&atom_3);
                           clipper::Coord_orth a4 = coot::co_gemmi(&atom_4);

                           double omega = clipper::Coord_orth::torsion(a1, a2, a3, a4);
                           double dist = clipper::Coord_orth::length(a2, a3);

                           bool add_it = true;
                           if (dist > 2.0) {
                              add_it = false;
                           }
                           if ((omega > 0.5 * M_PI) || (omega < -0.5 * M_PI)) {
                              // trans - ok
                           } else {
                              add_it = false; // cis
                           }
                           if (add_even_if_cis) {
                              add_it = true;
                           }

                           if (add_it) {
                              std::vector<bool> other_fixed_flags = make_fixed_flags(index1, index2, index3, index4);
                              for (std::size_t ii = 0; ii < other_fixed_flags.size(); ii++) {
                                 if (other_fixed_flags[ii])
                                    fixed_flags[ii] = true;
                              }

                              double target_omega = 180.0;
                              double esd = 2.0;
                              add(TRANS_PEPTIDE_RESTRAINT, index1, index2, index3, index4,
                                  fixed_flags,
                                  target_omega,
                                  esd,
                                  1.2, // dummy value
                                  1);
                              n_trans_peptide_torsion++;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   return n_trans_peptide_torsion;
}