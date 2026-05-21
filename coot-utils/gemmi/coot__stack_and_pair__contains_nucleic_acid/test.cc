#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <set>
#include "function.hh"

TEST(OracleTest, contains_nucleic_acid) {
   // Case: protein-only PDB
   {
      gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
      
      // Count atoms
      size_t atom_count = 0;
      for (const auto& model : st.models)
          for (const auto& chain : model.chains)
              for (const auto& res : chain.residues)
                  atom_count += res.atoms.size();
      
      EXPECT_EQ(atom_count, 2107);
      
      // Collect unique residues that have atoms
      std::set<const gemmi::Residue*> res_set;
      for (const auto& model : st.models)
          for (const auto& chain : model.chains)
              for (const auto& res : chain.residues)
                  res_set.insert(&res);
      
      // Check if any residue is a nucleotide
      bool has_nucleic_acid = false;
      for (const auto* res : res_set) {
          if (coot::util::is_nucleotide_gemmi(*res)) {
              has_nucleic_acid = true;
              break;
          }
      }
      EXPECT_FALSE(has_nucleic_acid);
   }

   // Case: empty selection
   {
      std::vector<gemmi::Residue> empty_residues;
      EXPECT_FALSE(coot::stack_and_pair::contains_nucleic_acid_gemmi(empty_residues));
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}