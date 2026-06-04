#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper to find a gemmi::Residue* by chain and residue sequence number
static gemmi::Residue *find_residue(gemmi::Structure &st, const std::string &chain_id, int seqnum) {
   if (st.models.empty()) return nullptr;
   for (gemmi::Chain &ch : st.models[0].chains) {
      if (ch.name != chain_id) continue;
      for (gemmi::Residue &res : ch.residues) {
         if (res.seqid.num.value == seqnum) return &res;
      }
   }
   return nullptr;
}

TEST(OracleTest, get_strand_sense) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: Two strands with >= 2 residues each (both forward in sequence)
    {
        std::vector<gemmi::Residue *> strand_1;
        std::vector<gemmi::Residue *> strand_2;

        gemmi::Residue *r10 = find_residue(st, "A", 10);
        gemmi::Residue *r11 = find_residue(st, "A", 11);
        gemmi::Residue *r20 = find_residue(st, "A", 20);
        gemmi::Residue *r21 = find_residue(st, "A", 21);

        if (r10) strand_1.push_back(r10);
        if (r11) strand_1.push_back(r11);
        if (r20) strand_2.push_back(r20);
        if (r21) strand_2.push_back(r21);

        auto sense = coot::secondary_structure_header_records::strand_relation_t::get_strand_sense_gemmi(strand_1, strand_2);

        EXPECT_EQ(strand_1.size(), 2u);
        EXPECT_EQ(strand_2.size(), 2u);
        EXPECT_EQ(static_cast<int>(sense), 2); // ANTI_PARALLEL
    }

    // Case 2: One strand with only 1 residue — should return NO_RESULT
    {
        std::vector<gemmi::Residue *> strand_1;
        std::vector<gemmi::Residue *> strand_2;

        gemmi::Residue *r10 = find_residue(st, "A", 10);
        gemmi::Residue *r20 = find_residue(st, "A", 20);
        gemmi::Residue *r21 = find_residue(st, "A", 21);

        if (r10) strand_1.push_back(r10);
        if (r20) strand_2.push_back(r20);
        if (r21) strand_2.push_back(r21);

        auto sense = coot::secondary_structure_header_records::strand_relation_t::get_strand_sense_gemmi(strand_1, strand_2);

        EXPECT_EQ(strand_1.size(), 1u);
        EXPECT_EQ(strand_2.size(), 2u);
        EXPECT_EQ(static_cast<int>(sense), 3); // NO_RESULT
    }

    // Case 3: Empty strands — should return NO_RESULT
    {
        std::vector<gemmi::Residue *> strand_1;
        std::vector<gemmi::Residue *> strand_2;

        auto sense = coot::secondary_structure_header_records::strand_relation_t::get_strand_sense_gemmi(strand_1, strand_2);

        EXPECT_EQ(strand_1.size(), 0u);
        EXPECT_EQ(strand_2.size(), 0u);
        EXPECT_EQ(static_cast<int>(sense), 3); // NO_RESULT
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}