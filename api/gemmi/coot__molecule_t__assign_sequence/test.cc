#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(AssignSequenceTest, BasicMutation) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());
    gemmi::Chain& chain = model.chains[0];
    
    std::string new_sequence(chain.residues.size(), 'A');
    bool success = coot::molecule_t::assign_sequence_gemmi(model, chain.name, new_sequence);
    EXPECT_TRUE(success);
    
    for (const auto& res : chain.residues) {
        EXPECT_EQ(res.name, "ALA");
    }
}

TEST(AssignSequenceTest, LengthMismatch) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());
    gemmi::Chain& chain = model.chains[0];
    
    std::string short_seq = "ABC";
    if (short_seq.length() == chain.residues.size()) {
        short_seq = "AB";
    }
    bool success = coot::molecule_t::assign_sequence_gemmi(model, chain.name, short_seq);
    EXPECT_FALSE(success);
}

TEST(AssignSequenceTest, SingleLetterToThreeLetter) {
    EXPECT_EQ(coot::molecule_t::single_letter_to_3_letter_code('A'), "ALA");
    EXPECT_EQ(coot::molecule_t::single_letter_to_3_letter_code('C'), "CYS");
    EXPECT_EQ(coot::molecule_t::single_letter_to_3_letter_code('G'), "GLY");
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}