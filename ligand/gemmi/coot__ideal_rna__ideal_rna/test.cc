#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, IdealRNA) {
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   {
      std::string RNA_or_DNA = "RNA";
      std::string form = "A";
      short single_stranded_flag = 1;
      std::string sequence = "aa";
      
      gemmi::Structure result = coot::ideal_rna::ideal_rna(RNA_or_DNA, form, single_stranded_flag, sequence, &st);
      
      int n_models = result.models.size();
      EXPECT_EQ(n_models, 1);
      if (n_models > 0) {
         const gemmi::Model& model = result.models[0];
         int nchains = model.chains.size();
         EXPECT_GT(nchains, 0);
         if (nchains > 0) {
            const gemmi::Chain& chain = model.chains[0];
            int nres = chain.residues.size();
            EXPECT_GT(nres, 0);
         }
      }
   }

   {
      std::string RNA_or_DNA = "DNA";
      std::string form = "B";
      short single_stranded_flag = 1;
      std::string sequence = "tt";
      
      gemmi::Structure result = coot::ideal_rna::ideal_rna(RNA_or_DNA, form, single_stranded_flag, sequence, &st);
      
      int n_models = result.models.size();
      EXPECT_EQ(n_models, 1);
      if (n_models > 0) {
         const gemmi::Model& model = result.models[0];
         int nchains = model.chains.size();
         EXPECT_GT(nchains, 0);
         if (nchains > 0) {
            const gemmi::Chain& chain = model.chains[0];
            int nres = chain.residues.size();
            EXPECT_GT(nres, 0);
         }
      }
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}