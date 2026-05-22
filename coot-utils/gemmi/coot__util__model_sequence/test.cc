#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(ModelSequenceTest, ModelSequence) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "load failed";

    // Get chain A
    const gemmi::Model& model = st.models[0];
    auto it = std::find_if(model.chains.begin(), model.chains.end(),
                           [](const gemmi::Chain& c) { return c.name == "A"; });
    ASSERT_NE(it, model.chains.end()) << "chain A not found";

    // Create the vector of pairs for model_sequence
    std::vector<std::pair<gemmi::Residue*, int>> sa;
    int idx = 0;
    for (const gemmi::Residue& res : it->residues) {
        sa.push_back(std::make_pair(const_cast<gemmi::Residue*>(&res), idx++));
    }

    // Test with allow_ligands = false
    {
        std::string seq = coot::util::model_sequence_gemmi(sa, false);
        EXPECT_EQ(sa.size(), 267u);
        EXPECT_EQ(seq.length(), 267u);
        EXPECT_EQ(seq, "MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRSTAIREISLLKELNHPNIVKLLDVIHTLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTEGAIKLADFLWFRAAEILLGCKYYSTAVDIWSLGCIYAEMVTRRALFPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLSQMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL");
    }

    // Test with allow_ligands = true
    {
        std::string seq = coot::util::model_sequence_gemmi(sa, true);
        EXPECT_EQ(sa.size(), 267u);
        EXPECT_EQ(seq.length(), 267u);
        EXPECT_EQ(seq, "MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRSTAIREISLLKELNHPNIVKLLDVIHTLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTEGAIKLADFLWFRAAEILLGCKYYSTAVDIWSLGCIYAEMVTRRALFPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLSQMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL");
    }

    // Edge case: empty vector
    {
        std::vector<std::pair<gemmi::Residue*, int>> empty_sa;
        std::string seq = coot::util::model_sequence_gemmi(empty_sa, false);
        EXPECT_EQ(seq.length(), 0u);
        EXPECT_EQ(seq, "");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}