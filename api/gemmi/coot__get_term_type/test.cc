#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_term_type) {
    // Load the example PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Helper lambda to find residue by chain ID and sequence number
    auto find_residue = [&](const std::string& chain_name, int seq_num) -> gemmi::CRA {
        for (auto& chain : model.chains) {
            if (chain.name == chain_name) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == seq_num) {
                        return gemmi::CRA{&chain, &res, nullptr};
                    }
                }
            }
        }
        return gemmi::CRA{nullptr, nullptr, nullptr};
    };

    // case: internal residue A/50 (has both neighbors) -> expected: "not-terminal-residue"
    {
        gemmi::CRA cra = find_residue("A", 50);
        ASSERT_NE(cra.residue, nullptr);
        std::string term_type = coot::get_term_type_gemmi(cra);
        EXPECT_EQ(term_type, "not-terminal-residue");
    }

    // case: C-terminus residue A/298 (no upstream neighbor) -> expected: "C"
    {
        gemmi::CRA cra = find_residue("A", 298);
        ASSERT_NE(cra.residue, nullptr);
        std::string term_type = coot::get_term_type_gemmi(cra);
        EXPECT_EQ(term_type, "C");
    }

    // case: N-terminus residue A/1 (no downstream neighbor) -> expected: "N"
    {
        gemmi::CRA cra = find_residue("A", 1);
        ASSERT_NE(cra.residue, nullptr);
        std::string term_type = coot::get_term_type_gemmi(cra);
        EXPECT_EQ(term_type, "N");
    }

    // case: singleton check for residue A/1 -> expected: "N"
    {
        gemmi::CRA cra = find_residue("A", 1);
        ASSERT_NE(cra.residue, nullptr);
        // Verify chain has 267 residues (matching original test)
        int nres = static_cast<int>(cra.chain->residues.size());
        EXPECT_EQ(nres, 267);
        std::string term_type = coot::get_term_type_gemmi(cra);
        EXPECT_EQ(term_type, "N");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}