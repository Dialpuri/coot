#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_atom_index_quad) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Get residue "A/10" (ALA in this structure)
    gemmi::Chain *chain = nullptr;
    for (auto &ch : st.models[0].chains) {
        if (ch.name == "A") {
            chain = &ch;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    gemmi::Residue *res = nullptr;
    for (auto &r : chain->residues) {
        if (r.seqid.num.value == 10) {
            res = &r;
            break;
        }
    }
    ASSERT_NE(res, nullptr);

    const std::string altconf("");

    // Case 1: torsion with atoms likely present (N-CA-C-CB)
    {
        coot::dict_torsion_restraint_t tr1("t1", "N", "CA", "C", "CB", 0.0, 0.0, 1);
        auto result1 = coot::get_atom_index_quad_gemmi(tr1, *res, altconf);

        EXPECT_TRUE(result1.first);
        EXPECT_EQ(result1.second.index1, 0);
        EXPECT_EQ(result1.second.index2, 1);
        EXPECT_EQ(result1.second.index3, 2);
        EXPECT_EQ(result1.second.index4, 4);
    }

    // Case 2: torsion with non-existent atom "OX1" — should fail
    {
        coot::dict_torsion_restraint_t tr2("t2", "N", "CA", "C", "OX1", 0.0, 0.0, 1);
        auto result2 = coot::get_atom_index_quad_gemmi(tr2, *res, altconf);

        EXPECT_FALSE(result2.first);
        EXPECT_EQ(result2.second.index1, 0);
        EXPECT_EQ(result2.second.index2, 1);
        EXPECT_EQ(result2.second.index3, 2);
        EXPECT_EQ(result2.second.index4, -1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}