#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, AddMainChainH) {
    // Read the test PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB file";

    gemmi::Model& model = st.models[0];
    ASSERT_GE(model.chains.size(), 1) << "No chains found";

    gemmi::Chain& chain = model.chains[0];

    // Find residues 10 and 11
    gemmi::Residue* res10 = nullptr;
    gemmi::Residue* res11 = nullptr;
    for (auto& res : chain.residues) {
        if (res.seqid.num.value == 10) res10 = &res;
        if (res.seqid.num.value == 11) res11 = &res;
    }

    ASSERT_NE(res10, nullptr) << "Failed to get residue 10";
    ASSERT_NE(res11, nullptr) << "Failed to get residue 11";

    // Case 1: residue 11 with previous residue 10
    {
        size_t h_count_before = 0;
        for (const auto& atom : res11->atoms) {
            if (atom.element.is_hydrogen()) h_count_before++;
        }

        coot::reduce::add_main_chain_H_gemmi(*res11, res10);

        size_t h_count_after = 0;
        for (const auto& atom : res11->atoms) {
            if (atom.element.is_hydrogen()) h_count_after++;
        }

        EXPECT_GT(h_count_after, 0) << "Residue 11 should have H after adding";
    }

    // Case 2: residue 10 with nullptr as previous residue
    {
        size_t h_count_before = 0;
        for (const auto& atom : res10->atoms) {
            if (atom.element.is_hydrogen()) h_count_before++;
        }

        coot::reduce::add_main_chain_H_gemmi(*res10, nullptr);

        size_t h_count_after = 0;
        for (const auto& atom : res10->atoms) {
            if (atom.element.is_hydrogen()) h_count_after++;
        }

        EXPECT_EQ(h_count_after, h_count_before) << "Residue 10 should have same H count with nullptr";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}