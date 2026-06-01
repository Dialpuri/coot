#include <gtest/gtest.h>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, residues_sort_function) {
    // Setup: load CIF and get residues
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    gemmi::Model& model = st.models[0];

    gemmi::CRA cra1{};
    gemmi::CRA cra2{};

    for (auto& chain : model.chains) {
        if (chain.name == "A") {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 1) {
                    cra1 = gemmi::CRA{&chain, &res, nullptr};
                }
                if (res.seqid.num.value == 10) {
                    cra2 = gemmi::CRA{&chain, &res, nullptr};
                }
            }
        }
    }

    EXPECT_NE(cra1.residue, nullptr);
    EXPECT_NE(cra2.residue, nullptr);

    // Case 1: Normal comparison - r1 (seq 1) should come before r2 (seq 10)
    {
        bool result = coot::util::residues_sort_function_gemmi(cra1, cra2);
        EXPECT_TRUE(result);
    }

    // Case 2: Reversed order - r2 (seq 10) should NOT come before r1 (seq 1)
    {
        bool result = coot::util::residues_sort_function_gemmi(cra2, cra1);
        EXPECT_FALSE(result);
    }

    // Case 3: Both nullptr
    {
        bool result = coot::util::residues_sort_function_gemmi({}, {});
        EXPECT_FALSE(result);
    }

    // Case 4: r1 is nullptr, r2 is valid
    {
        bool result = coot::util::residues_sort_function_gemmi({}, cra2);
        EXPECT_TRUE(result);
    }

    // Case 5: r1 is valid, r2 is nullptr
    {
        bool result = coot::util::residues_sort_function_gemmi(cra1, {});
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}