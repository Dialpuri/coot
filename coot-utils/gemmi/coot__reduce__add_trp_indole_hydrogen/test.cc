#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(GemmiReduceTest, add_trp_indole_hydrogen) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model& model = st.models[0];

    // case 1: LEU 166 has backbone atoms N, CA, C — should add H
    {
        gemmi::Residue* res = nullptr;
        for (auto& chain : model.chains) {
            if (chain.name != "A") continue;
            for (auto& residue : chain.residues) {
                if (residue.seqid.num.value == 166) {
                    res = &residue;
                    break;
                }
            }
        }
        ASSERT_NE(res, nullptr);

        const gemmi::Chain* chain_a = model.find_chain("A");
        ASSERT_NE(chain_a, nullptr);
        EXPECT_STREQ(chain_a->name.c_str(), "A");
        EXPECT_EQ(res->seqid.num.value, 166);

        int atoms_before = static_cast<int>(res->atoms.size());
        EXPECT_EQ(atoms_before, 21);

        coot::reduce::add_trp_indole_hydrogen_gemmi("HN1", " N  ", " CA ", " C  ", 1.01, *res);

        int atoms_after = static_cast<int>(res->atoms.size());
        EXPECT_EQ(atoms_after, 22);
    }

    // case 2: TRP 167 missing side chain atoms — should add nothing
    {
        gemmi::Residue* res = nullptr;
        for (auto& chain : model.chains) {
            if (chain.name != "A") continue;
            for (auto& residue : chain.residues) {
                if (residue.seqid.num.value == 167) {
                    res = &residue;
                    break;
                }
            }
        }
        ASSERT_NE(res, nullptr);

        const gemmi::Chain* chain_a_2 = model.find_chain("A");
        ASSERT_NE(chain_a_2, nullptr);
        EXPECT_STREQ(chain_a_2->name.c_str(), "A");
        EXPECT_EQ(res->seqid.num.value, 167);

        int atoms_before = static_cast<int>(res->atoms.size());
        EXPECT_EQ(atoms_before, 9);

        coot::reduce::add_trp_indole_hydrogen_gemmi("HD1", " CG ", " CD1", " NE1", 1.01, *res);

        int atoms_after = static_cast<int>(res->atoms.size());
        EXPECT_EQ(atoms_after, 9);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}