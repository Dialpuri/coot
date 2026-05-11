#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, density_box_t_default_construction) {
   // case: default construction
   coot::density_box_t db;
   EXPECT_TRUE(db.empty());
   EXPECT_EQ(db.density_box, nullptr);
   EXPECT_EQ(db.residue_p, nullptr);
   EXPECT_EQ(db.n_steps, 0);
   EXPECT_NEAR(db.mean, 0.0, 1e-4);
   EXPECT_NEAR(db.var, -1.0, 1e-4);
   EXPECT_FALSE(db.is_weird);
}

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}