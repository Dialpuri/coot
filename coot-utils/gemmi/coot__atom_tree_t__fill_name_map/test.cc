#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atom_tree_t_fill_name_map) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get first model
    EXPECT_GE(st.models.size(), 1);
    gemmi::Model* model = &st.models[0];

    // Get first chain, first residue
    EXPECT_GE(model->chains.size(), 1);
    gemmi::Chain* chain = &model->chains[0];
    EXPECT_GE(chain->residues.size(), 1);

    // Case 1: Use first residue (MET 1)
    gemmi::Residue* res = &chain->residues[0];
    ASSERT_NE(res, nullptr);

    EXPECT_EQ(chain->name, "A");
    EXPECT_EQ(res->seqid.num.value, 1);
    EXPECT_EQ(res->name, "MET");

    // Case 1: Use first residue with empty altconf
    std::string altconf1 = "";

    auto name_to_index = coot::fill_name_map_gemmi(*res, altconf1);

    int size_after = static_cast<int>(name_to_index.size());
    EXPECT_EQ(size_after, 8);

    // Verify specific entries from name_to_index map (oracle ground truth)
    // Keys have leading/trailing spaces as per PDB atom name format
    EXPECT_EQ(name_to_index.at(" CA ").index(), 1);
    EXPECT_EQ(name_to_index.at(" C  ").index(), 2);
    EXPECT_EQ(name_to_index.at(" CB ").index(), 4);
    EXPECT_EQ(name_to_index.at(" CG ").index(), 5);
    EXPECT_EQ(name_to_index.at(" CE ").index(), 7);

    // Case 2: non-matching altconf "X" (should map only atoms with altLoc=="X", likely none)
    std::string altconf2 = "X";

    auto name_to_index2 = coot::fill_name_map_gemmi(*res, altconf2);

    int size_after2 = static_cast<int>(name_to_index2.size());
    EXPECT_EQ(size_after2, 8);

    // Verify same entries as Case 1
    EXPECT_EQ(name_to_index2.at(" CA ").index(), 1);
    EXPECT_EQ(name_to_index2.at(" C  ").index(), 2);
    EXPECT_EQ(name_to_index2.at(" CB ").index(), 4);
    EXPECT_EQ(name_to_index2.at(" CG ").index(), 5);
    EXPECT_EQ(name_to_index2.at(" CE ").index(), 7);

    // Case 3: Now manually call fill_name_map on same object with different altconf
    // Clear the map first, then call with empty altconf again
    name_to_index2.clear();
    std::string altconf3 = "";

    int size_before = static_cast<int>(name_to_index2.size());
    EXPECT_EQ(size_before, 0);

    name_to_index2 = coot::fill_name_map_gemmi(*res, altconf3);

    int size_after3 = static_cast<int>(name_to_index2.size());
    EXPECT_EQ(size_after3, 8);

    // Verify entries after fill_name_map
    EXPECT_EQ(name_to_index2.at(" CA ").index(), 1);
    EXPECT_EQ(name_to_index2.at(" C  ").index(), 2);
    EXPECT_EQ(name_to_index2.at(" CB ").index(), 4);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}