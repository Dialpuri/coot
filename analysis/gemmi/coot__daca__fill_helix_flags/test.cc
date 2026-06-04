#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, fill_helix_flags) {
    // Load structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    int n_models = static_cast<int>(st.models.size());
    EXPECT_EQ(n_models, 1);

    if (n_models < 1) {
        FAIL() << "No models loaded";
    }

    gemmi::Model *model = &st.models[0];
    ASSERT_NE(model, nullptr);

    EXPECT_EQ(model->chains.size(), 2);

    // Case 1: valid model — call fill_helix_flags
    {
        std::vector<gemmi::CRA> helical_residues;

        // Clear helical_residues before call (ensure clean state)
        helical_residues.clear();
        size_t before = helical_residues.size();

        // Call the function
        coot::daca::fill_helix_flags_gemmi(model, helical_residues);

        size_t after = helical_residues.size();

        EXPECT_EQ(before, 0u);
        EXPECT_EQ(after, 54u);

        // Print some details about the helical residues found
        if (after > 0) {
            std::string sample_residues;
            for (unsigned int i = 0; i < after && i < 10; i++) {
                if (i > 0) sample_residues += ", ";
                gemmi::CRA cra = helical_residues[i];
                ASSERT_NE(cra.residue, nullptr);
                ASSERT_NE(cra.chain, nullptr);
                sample_residues += cra.chain->name + "/" +
                                  std::to_string(cra.residue->seqid.num.value) + "/" +
                                  cra.residue->name;
            }
            EXPECT_EQ(sample_residues, "A/47/THR, A/48/ALA, A/49/ILE, A/50/ARG, A/51/GLU, A/52/ILE, A/53/SER, A/54/LEU, A/88/LYS, A/89/LYS");
        }
    }

    // Case 2: null model — should be a no-op
    {
        std::vector<gemmi::CRA> helical_residues;
        helical_residues.clear();
        size_t before = helical_residues.size();

        coot::daca::fill_helix_flags_gemmi(nullptr, helical_residues);

        size_t after = helical_residues.size();

        EXPECT_EQ(before, 0u);
        EXPECT_EQ(after, 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}