#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, deep_copy_this_residue_with_atom_index_and_afix_transfer) {
    // Load source molecule
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get a residue from chain A, seqnum 1 (MET)
    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());

    // Find the first residue (chain A, MET, seqnum 1)
    gemmi::Chain* chain_ptr = &model.chains[0];
    EXPECT_EQ(chain_ptr->name, "A");

    ASSERT_FALSE(chain_ptr->residues.empty());
    gemmi::Residue* source_residue = &chain_ptr->residues[0];
    EXPECT_EQ(source_residue->seqid.num.value, 1);
    EXPECT_EQ(source_residue->name, "MET");

    int source_atom_count = static_cast<int>(source_residue->atoms.size());
    EXPECT_EQ(source_atom_count, 8);

    // === Case 1: whole_residue_flag=1 (copy all atoms) ===
    {
        std::string altconf("");
        short int whole_residue_flag = 1;

        gemmi::CRA cra{chain_ptr, source_residue, nullptr};
        auto result = coot::util::deep_copy_this_residue_with_atom_index_and_afix_transfer_gemmi(
            cra, altconf, whole_residue_flag);

        int result_atom_count = result ? static_cast<int>(result->atoms.size()) : 0;

        EXPECT_NE(result, nullptr);  // residue_created: true
        EXPECT_EQ(result_atom_count, 8);   // copied_atoms: 8
        EXPECT_EQ(result ? result->seqid.num.value : -1, 1);  // result_seqnum: 1
    }

    // === Case 2: whole_residue_flag=0, altconf="" (copy only atoms with altLoc "" ) ===
    {
        std::string altconf("");
        short int whole_residue_flag = 0;

        gemmi::CRA cra{chain_ptr, source_residue, nullptr};
        auto result = coot::util::deep_copy_this_residue_with_atom_index_and_afix_transfer_gemmi(
            cra, altconf, whole_residue_flag);

        int result_atom_count = result ? static_cast<int>(result->atoms.size()) : 0;

        EXPECT_NE(result, nullptr);  // residue_created: true
        EXPECT_EQ(result_atom_count, 8);   // copied_atoms: 8
        EXPECT_EQ(result ? result->seqid.num.value : -1, 1);  // result_seqnum: 1
    }

    // === Case 3: whole_residue_flag=0, altconf="X" (non-matching altLoc -> fewer atoms) ===
    {
        std::string altconf("X");
        short int whole_residue_flag = 0;

        gemmi::CRA cra{chain_ptr, source_residue, nullptr};
        auto result = coot::util::deep_copy_this_residue_with_atom_index_and_afix_transfer_gemmi(
            cra, altconf, whole_residue_flag);

        int result_atom_count = result ? static_cast<int>(result->atoms.size()) : 0;

        EXPECT_NE(result, nullptr);  // residue_created: true
        EXPECT_EQ(result_atom_count, 8);   // copied_atoms: 8
        EXPECT_EQ(result ? result->seqid.num.value : -1, 1);  // result_seqnum: 1
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}