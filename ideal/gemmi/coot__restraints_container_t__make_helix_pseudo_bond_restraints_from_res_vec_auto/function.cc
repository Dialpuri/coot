#include "function.hh"

#include <algorithm>
#include <vector>
#include <string>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__compare_to_helix__188839/gemmi/function.hh"

namespace coot {

int make_helix_pseudo_bond_restraints_from_res_vec_auto_gemmi(
    const std::vector<std::pair<bool, gemmi::Residue*>>& residues_vec
) {
   float pseudo_bond_esd = 0.1;
   unsigned int n_helical_restraints = 0;

   // Collect non-null residues
   std::vector<gemmi::Residue*> sorted_residues;
   sorted_residues.reserve(residues_vec.size());
   for (unsigned int i = 0; i < residues_vec.size(); i++)
      if (residues_vec[i].second)
         sorted_residues.push_back(residues_vec[i].second);

   // Sort by chain name then sequence number
   std::sort(sorted_residues.begin(), sorted_residues.end(),
       [](const gemmi::Residue* a, const gemmi::Residue* b) {
           // gemmi::Residue doesn't have a chain pointer, so sort by seqid only
           return a->seqid.num.value < b->seqid.num.value;
       });

   for (unsigned int i = 0; i < sorted_residues.size(); i++) {

      if ((i + 3) >= sorted_residues.size()) continue;

      // Build test_helical_residues with up to 4 consecutive residues
      std::vector<gemmi::Residue*> test_helical_residues;
      test_helical_residues.reserve(5);

      gemmi::Residue* residue_0 = sorted_residues[i];
      for (unsigned int iir = 0; iir < 4; iir++) {
         if ((i + iir) < sorted_residues.size()) {
            test_helical_residues.push_back(sorted_residues[i + iir]);
         }
      }

      // Check sane residue numbers
      bool sane_residue_numbers = false;
      if (test_helical_residues.size() == 4)
         if ((test_helical_residues[0]->seqid.num.value + 3) == test_helical_residues[3]->seqid.num.value)
            sane_residue_numbers = true;

      if (test_helical_residues.size() >= 5)
         if ((test_helical_residues[0]->seqid.num.value + 4) == test_helical_residues[4]->seqid.num.value)
            sane_residue_numbers = true;

      // Call the gemmi port of compare_to_helix (single-arg overload that uses default ref positions)
      helical_results_t hr = coot::compare_to_helix_gemmi(test_helical_residues);

      if (hr.is_alpha_helix_like && sane_residue_numbers) {

         int index_1 = -1; // O atom index
         int index_2 = -1; // N (n+4) atom index
         int index_3 = -1; // N (n+3) atom index
         const gemmi::Atom* at_1 = nullptr;
         const gemmi::Atom* at_2 = nullptr;
         const gemmi::Atom* at_3 = nullptr;

         bool do_i_plus_4_also = true;
         bool do_i_plus_3 = true;

         if ((i + 4) >= sorted_residues.size())
            do_i_plus_4_also = false;

         if (do_i_plus_4_also) {
            if (sorted_residues[i + 4]->seqid.num.value != (sorted_residues[i]->seqid.num.value + 4))
               do_i_plus_4_also = false;
         }

         if ((i + 3) < sorted_residues.size()) {
            if (sorted_residues[i + 3]->seqid.num.value != (sorted_residues[i]->seqid.num.value + 3))
               do_i_plus_3 = false;
         }

         // Find O atom in residue[i] (by element name "O")
         for (const gemmi::Atom& atom : sorted_residues[i]->atoms) {
            if (atom.element.name() == "O") {
               at_1 = &atom;
               break;
            }
         }

         // Find N atom in residue[i+3]
         if (do_i_plus_3 && (i + 3) < sorted_residues.size()) {
            for (const gemmi::Atom& atom : sorted_residues[i + 3]->atoms) {
               if (atom.element.name() == "N") {
                  at_3 = &atom;
                  break;
               }
            }
         }

         // Find N atom in residue[i+4]
         if (do_i_plus_4_also && (i + 4) < sorted_residues.size()) {
            for (const gemmi::Atom& atom : sorted_residues[i + 4]->atoms) {
               if (atom.element.name() == "N") {
                  at_2 = &atom;
                  break;
               }
            }
         }

         if (at_1 && at_2 && at_3) {
            // Would add two restraints in the original (i+3 and i+4)
            n_helical_restraints += 2;
         } else if (at_1 && at_3) {
            // Would add one restraint (i+3 only)
            n_helical_restraints += 1;
         }
      }
   }

   return static_cast<int>(n_helical_restraints);
}

} // namespace coot