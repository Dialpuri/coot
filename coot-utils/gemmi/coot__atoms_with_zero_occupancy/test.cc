#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atoms_with_zero_occupancy) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: Normal PDB — all atoms have occupancy ~1.0, so expect 0 zero-occupancy atoms
    {
        std::vector<gemmi::Atom*> result = coot::atoms_with_zero_occupancy_gemmi(st.models[0]);
        EXPECT_EQ(result.size(), 0u);
    }

    // Case 2: Set some atoms to zero occupancy, then call the function
    {
        gemmi::Model& model = st.models[0];
        ASSERT_FALSE(model.chains.empty());

        gemmi::Chain& chain = model.chains[0];
        ASSERT_FALSE(chain.residues.empty());

        gemmi::Residue& res = chain.residues[0];
        ASSERT_FALSE(res.atoms.empty());

        gemmi::Atom* target_atom = nullptr;
        double saved_occ = 0;

        for (auto& at : res.atoms) {
            if (at.occ >= 0.01) {
                saved_occ = at.occ;
                at.occ = 0.0;
                target_atom = &at;
                break;
            }
        }
        ASSERT_NE(target_atom, nullptr);

        {
            std::vector<gemmi::Atom*> result = coot::atoms_with_zero_occupancy_gemmi(st.models[0]);
            EXPECT_EQ(result.size(), 1u);
            if (result.size() > 0) {
                EXPECT_EQ(result[0], target_atom);
            }
        }

        // Restore
        target_atom->occ = saved_occ;
    }

    // Case 3: Set multiple atoms to zero occupancy
    {
        gemmi::Model& model = st.models[0];
        ASSERT_FALSE(model.chains.empty());

        int zero_count = 0;
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                for (auto& at : res.atoms) {
                    if (at.occ >= 0.01) {
                        at.occ = 0.0;
                        zero_count++;
                    }
                }
            }
        }

        std::vector<gemmi::Atom*> result = coot::atoms_with_zero_occupancy_gemmi(st.models[0]);
        EXPECT_EQ(result.size(), 2107u);
        EXPECT_EQ(result.size(), static_cast<size_t>(zero_count));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}