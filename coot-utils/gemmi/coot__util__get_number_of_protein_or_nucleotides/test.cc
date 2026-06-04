#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_number_of_protein_or_nucleotides) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_EQ(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];

    // Case 1: chain A — protein residues (transferase 2VTQ)
    {
        gemmi::Chain* chain_p = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                chain_p = &chain;
                break;
            }
        }
        EXPECT_NE(chain_p, nullptr);

        std::pair<unsigned int, unsigned int> result =
            coot::util::get_number_of_protein_or_nucleotides_gemmi(chain_p);

        EXPECT_EQ(result.first, 267u);   // protein_count
        EXPECT_EQ(result.second, 0u);    // nucleotide_count
    }

    // Case 2: chain B — only HOH water, 0 protein/0 nucleotide
    {
        gemmi::Chain* chain_p = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "B") {
                chain_p = &chain;
                break;
            }
        }
        EXPECT_NE(chain_p, nullptr);

        std::pair<unsigned int, unsigned int> result =
            coot::util::get_number_of_protein_or_nucleotides_gemmi(chain_p);

        EXPECT_EQ(result.first, 0u);     // protein_count
        EXPECT_EQ(result.second, 0u);    // nucleotide_count
    }

    // Case 3: null chain pointer — tests the guard clause
    {
        gemmi::Chain* chain_p = nullptr;

        std::pair<unsigned int, unsigned int> result =
            coot::util::get_number_of_protein_or_nucleotides_gemmi(chain_p);

        EXPECT_EQ(result.first, 0u);     // protein_count
        EXPECT_EQ(result.second, 0u);    // nucleotide_count
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}