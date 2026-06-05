#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, restraints_container_add_link_angle) {
    coot::protein_geometry geom;
    geom.init_standard();

    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

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

    // Case 1: valid peptide link "TRANS" between two consecutive residues
    {
        std::string link_type = "TRANS";
        short int is_fixed_first = 0;
        short int is_fixed_second = 0;

        int nangle = coot::add_link_angle_gemmi(link_type, cra_res1, cra_res2,
                                                is_fixed_first, is_fixed_second,
                                                geom);

        EXPECT_EQ(nangle, 4);
    }

    // Case 2: non-existent link type — should return 0
    {
        std::string link_type = "NONEXISTENT_LINK";
        short int is_fixed_first = 0;
        short int is_fixed_second = 0;

        int nangle = coot::add_link_angle_gemmi(link_type, cra_res1, cra_res2,
                                                is_fixed_first, is_fixed_second,
                                                geom);

        EXPECT_EQ(nangle, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}