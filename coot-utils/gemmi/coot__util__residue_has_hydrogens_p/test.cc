#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ResidueHasHydrogens) {
    // case 1: residue with hydrogens (MET 1 has H atoms)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
        gemmi::Model& model = st.models[0];
        
        // Find chain A, residue 1 (MET)
        const gemmi::Residue* res1 = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 1) {
                        res1 = &res;
                        break;
                    }
                }
            }
        }
        
        EXPECT_TRUE(res1 != nullptr) << "residue should be found";
        if (res1) {
            EXPECT_EQ(res1->name, "MET");
            EXPECT_EQ(1, res1->seqid.num.value);
            // chain.name is already verified by the outer loop condition
        }
        EXPECT_TRUE(coot::util::residue_has_hydrogens_p_gemmi(*res1)) << "MET residue should have hydrogens";
    }

    // case 2: another residue with hydrogens (GLU 2)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
        gemmi::Model& model = st.models[0];
        
        // Find chain A, residue 2 (GLU)
        const gemmi::Residue* res2 = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 2) {
                        res2 = &res;
                        break;
                    }
                }
            }
        }
        
        EXPECT_TRUE(res2 != nullptr) << "residue should be found";
        if (res2) {
            EXPECT_EQ(res2->name, "GLU");
            EXPECT_EQ(2, res2->seqid.num.value);
        }
        EXPECT_TRUE(coot::util::residue_has_hydrogens_p_gemmi(*res2)) << "GLU residue should have hydrogens";
    }

    // case 3: null residue (non-existent residue) — in gemmi we cannot have a null pointer,
    // so this test case is not directly translatable. We skip this case as it has no equivalent in gemmi's design.
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}