#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_atom_index_pairs) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get residue //A/1 (MET)
    const gemmi::Residue *res = nullptr;
    std::string chain_id;
    for (gemmi::Model &model : st.models) {
        for (gemmi::Chain &chain : model.chains) {
            for (gemmi::Residue &r : chain.residues) {
                if (r.seqid.num.value == 1 && chain.name == "A") {
                    res = &r;
                    chain_id = chain.name;
                    break;
                }
            }
        }
    }
    ASSERT_NE(res, nullptr);
    int seqnum = res->seqid.num.value;
    std::string resname = res->name;
    int n_atoms = static_cast<int>(res->atoms.size());

    EXPECT_EQ(chain_id, "A");
    EXPECT_EQ(seqnum, 1);
    EXPECT_EQ(resname, "MET");
    EXPECT_EQ(n_atoms, 8);

    // === Oracle cases (all use unpadded names → early-return guard, size 0) ===

    // Case 1: valid atom pairs in MET (N-CA and CA-C) — oracle outputs size 0
    {
        std::vector<coot::atom_name_pair> pairs;
        pairs.push_back(coot::atom_name_pair("N", "CA"));
        pairs.push_back(coot::atom_name_pair("CA", "C"));

        std::vector<coot::atom_index_pair> result = coot::monomer_utils::get_atom_index_pairs_gemmi(pairs, *res);

        EXPECT_EQ(result.size(), 0u);
    }

    // Case 2: single pair with atoms that exist in MET — oracle outputs size 0
    {
        std::vector<coot::atom_name_pair> pairs;
        pairs.push_back(coot::atom_name_pair("CB", "SG"));

        std::vector<coot::atom_index_pair> result = coot::monomer_utils::get_atom_index_pairs_gemmi(pairs, *res);

        EXPECT_EQ(result.size(), 0u);
    }

    // Case 3: pair with first atom not in residue — oracle outputs size 0
    {
        std::vector<coot::atom_name_pair> pairs;
        pairs.push_back(coot::atom_name_pair("NONEXIST", "CA"));

        std::vector<coot::atom_index_pair> result = coot::monomer_utils::get_atom_index_pairs_gemmi(pairs, *res);

        EXPECT_EQ(result.size(), 0u);
    }

    // === Complementary cases (padded 4-char names exercise core logic path) ===
    // MET atom order: 0:" N  ", 1:" CA ", 2:" C  ", 3:" O  ", 4:" CB ", 5:" CG ", 6:" SD ", 7:" CE "

    // Case 4: two padded atom pairs (N-CA and CA-C)
    {
        std::vector<coot::atom_name_pair> pairs;
        pairs.push_back(coot::atom_name_pair(" N  ", " CA "));
        pairs.push_back(coot::atom_name_pair(" CA ", " C  "));

        std::vector<coot::atom_index_pair> result = coot::monomer_utils::get_atom_index_pairs_gemmi(pairs, *res);

        EXPECT_EQ(result.size(), 2u);
        EXPECT_EQ(result[0].index1, 0);  // N
        EXPECT_EQ(result[0].index2, 1);  // CA
        EXPECT_EQ(result[1].index1, 1);  // CA
        EXPECT_EQ(result[1].index2, 2);  // C
    }

    // Case 5: single padded pair CB-CG in MET
    {
        std::vector<coot::atom_name_pair> pairs;
        pairs.push_back(coot::atom_name_pair(" CB ", " CG "));

        std::vector<coot::atom_index_pair> result = coot::monomer_utils::get_atom_index_pairs_gemmi(pairs, *res);

        EXPECT_EQ(result.size(), 1u);
        EXPECT_EQ(result[0].index1, 4);  // CB
        EXPECT_EQ(result[0].index2, 5);  // CG
    }

    // Case 6: padded pair with non-existent second atom → returns 0
    {
        std::vector<coot::atom_name_pair> pairs;
        pairs.push_back(coot::atom_name_pair(" N  ", " NONEXIST "));

        std::vector<coot::atom_index_pair> result = coot::monomer_utils::get_atom_index_pairs_gemmi(pairs, *res);

        EXPECT_EQ(result.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}