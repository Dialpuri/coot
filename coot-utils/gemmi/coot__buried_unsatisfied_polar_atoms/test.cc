#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, buried_unsatisfied_polar_atoms) {
    // Case: standard protein with hydrogens (example-hydrogen.pdb)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
        EXPECT_GE(st.models.size(), 1u);

        int n_atoms_before = 0;
        for (gemmi::Chain& chain : st.models[0].chains) {
            for (gemmi::Residue& res : chain.residues) {
                n_atoms_before += static_cast<int>(res.atoms.size());
            }
        }

        EXPECT_EQ(n_atoms_before, 4280);

        // Call the function — it is void and prints "selected X atoms" internally
        EXPECT_NO_THROW(coot::buried_unsatisfied_polar_atoms_gemmi(st));

        int n_atoms_after = 0;
        for (gemmi::Chain& chain : st.models[0].chains) {
            for (gemmi::Residue& res : chain.residues) {
                n_atoms_after += static_cast<int>(res.atoms.size());
            }
        }

        EXPECT_EQ(n_atoms_after, 4280);
        EXPECT_TRUE(n_atoms_before == n_atoms_after);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}