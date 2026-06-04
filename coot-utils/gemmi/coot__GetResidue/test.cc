#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetResidue) {
   // Case 1: residue with multiple atoms
   {
      coot::minimol::residue res_in(10, "ALA");
      res_in.ins_code = "";
      res_in.addatom("N",   "N", 1.0f, 2.0f, 3.0f, "", 20.0f, 1.0f);
      res_in.addatom("CA",  "C", 2.0f, 3.0f, 4.0f, "", 22.0f, 1.0f);
      res_in.addatom("C",   "C", 3.0f, 4.0f, 5.0f, "", 24.0f, 1.0f);
      res_in.addatom("O",   "O", 4.0f, 5.0f, 6.0f, "", 26.0f, 1.0f);

      EXPECT_EQ(res_in.name, "ALA");
      EXPECT_EQ(res_in.seqnum, 10);
      EXPECT_EQ(res_in.ins_code, "");
      EXPECT_EQ(res_in.atoms.size(), 4u);

      gemmi::Residue res = coot::GetResidue_gemmi(res_in);

      EXPECT_EQ(res.name, "ALA");
      EXPECT_EQ(res.seqid.num, 10);
      EXPECT_EQ(res.seqid.icode, ' ');
      EXPECT_EQ(res.atoms.size(), 4u);

      // Print first atom details
      if (res.atoms.size() > 0) {
         const auto &at = res.atoms[0];
         EXPECT_EQ(at.name, "N");
         EXPECT_EQ(std::string(at.element.name()), "N");
         EXPECT_NEAR(at.pos.x, 1.0f, 1e-4);
         EXPECT_NEAR(at.pos.y, 2.0f, 1e-4);
         EXPECT_NEAR(at.pos.z, 3.0f, 1e-4);
      }
   }

   // Case 2: empty residue (no atoms)
   {
      coot::minimol::residue res_in(42, "GLY");
      res_in.ins_code = "A";

      EXPECT_EQ(res_in.name, "GLY");
      EXPECT_EQ(res_in.seqnum, 42);
      EXPECT_EQ(res_in.ins_code, "A");
      EXPECT_EQ(res_in.atoms.size(), 0u);

      gemmi::Residue res = coot::GetResidue_gemmi(res_in);

      EXPECT_EQ(res.name, "GLY");
      EXPECT_EQ(res.seqid.num, 42);
      EXPECT_EQ(res.seqid.icode, 'A');
      EXPECT_EQ(res.atoms.size(), 0u);
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}