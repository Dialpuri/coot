#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_link_torsion_gemmi) {
    coot::protein_geometry pg;
    pg.init_standard();

    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model &model = st.models[0];

    // Find residues 1 and 2 in chain A
    gemmi::CRA cra_res1{nullptr, nullptr, nullptr};
    gemmi::CRA cra_res2{nullptr, nullptr, nullptr};

    for (auto &chain : model.chains) {
        if (chain.name == "A") {
            for (auto &res : chain.residues) {
                if (res.seqid.num.value == 1 && cra_res1.residue == nullptr) {
                    cra_res1 = gemmi::CRA{&chain, &res, nullptr};
                }
                if (res.seqid.num.value == 2 && cra_res2.residue == nullptr) {
                    cra_res2 = gemmi::CRA{&chain, &res, nullptr};
                }
            }
        }
    }

    ASSERT_NE(cra_res1.residue, nullptr);
    ASSERT_NE(cra_res2.residue, nullptr);

    std::string link_type = "TRANS";

    // Case 1: valid adjacent residues
    {
        short int is_fixed_first = 0;
        short int is_fixed_second = 0;

        int result = coot::add_link_torsion_gemmi(link_type, cra_res1, cra_res2,
                                                   is_fixed_first, is_fixed_second, pg);

        EXPECT_EQ(result, 0);
    }

    // Case 2: same residue passed twice (edge case)
    {
        short int is_fixed_first = 0;
        short int is_fixed_second = 0;

        int result = coot::add_link_torsion_gemmi(link_type, cra_res1, cra_res1,
                                                   is_fixed_first, is_fixed_second, pg);

        EXPECT_EQ(result, 0);
    }

    // Case 3: non-existent link type
    {
        std::string bad_link_type = "NONEXISTENT";
        short int is_fixed_first = 0;
        short int is_fixed_second = 0;

        int result = coot::add_link_torsion_gemmi(bad_link_type, cra_res1, cra_res2,
                                                   is_fixed_first, is_fixed_second, pg);

        EXPECT_EQ(result, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}