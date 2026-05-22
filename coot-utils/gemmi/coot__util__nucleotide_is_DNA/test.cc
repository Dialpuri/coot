#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, nucleotide_is_DNA) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Find residue //A/10
    const gemmi::Residue* res = nullptr;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& r : chain.residues) {
                if (chain.name == "A" && r.seqid.num.value == 10) {
                    res = &r;
                    break;
                }
            }
            if (res) break;
        }
        if (res) break;
    }
    
    // Case: example.pdb has no nucleotides, so we may not find residue 10
    // or it might be a protein residue. According to the oracle, nullptr
    // or non-nucleotide should return true (DNA logic - no O2' found)
    if (!res) {
        bool result = coot::util::nucleotide_is_DNA_gemmi(nullptr);
        EXPECT_TRUE(result) << "nullptr residue should return true (DNA logic)";
        return;
    }
    
    // If a residue exists (edge case), test it
    {
        bool result = coot::util::nucleotide_is_DNA_gemmi(res);
        EXPECT_TRUE(result) << "Residue //A/10 should return true (DNA logic)";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}