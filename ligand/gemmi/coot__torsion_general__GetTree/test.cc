#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, torsion_general) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Case 1: Residue A/10 (ILE) - 8 atoms
    {
        const gemmi::Chain* chain = model.find_chain("A");
        ASSERT_NE(chain, nullptr);
        const gemmi::Residue* res = nullptr;
        for (const auto& r : chain->residues) {
            if (r.seqid.num.value == 10 && r.name == "ILE") { res = &r; break; }
        }
        ASSERT_NE(res, nullptr);

        EXPECT_EQ(res->seqid.num.value, 10);
        EXPECT_EQ(res->atoms.size(), 8u);

        Tree tree = coot::torsion_general::GetTree_gemmi(*res);
        EXPECT_EQ(tree.GetNumberOfVertices(), 8);
    }

    // Case 2: Residue A/20 (LYS) - 5 atoms
    {
        const gemmi::Chain* chain2 = model.find_chain("A");
        ASSERT_NE(chain2, nullptr);
        const gemmi::Residue* res2 = nullptr;
        for (const auto& r : chain2->residues) {
            if (r.seqid.num.value == 20 && r.name == "LYS") { res2 = &r; break; }
        }
        ASSERT_NE(res2, nullptr);

        EXPECT_EQ(res2->seqid.num.value, 20);
        EXPECT_EQ(res2->atoms.size(), 5u);

        Tree tree2 = coot::torsion_general::GetTree_gemmi(*res2);
        EXPECT_EQ(tree2.GetNumberOfVertices(), 5);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}