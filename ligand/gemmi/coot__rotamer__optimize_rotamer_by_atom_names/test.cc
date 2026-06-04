#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, optimize_rotamer_by_atom_names) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Find a suitable residue (PHE, ASP, GLU, or TYR)
    gemmi::CRA target_cra{nullptr, nullptr, nullptr};
    for (auto& chain : model.chains) {
        for (auto& res : chain.residues) {
            std::string rname = res.name;
            if (rname == "PHE" || rname == "ASP" || rname == "GLU" || rname == "TYR") {
                target_cra = gemmi::CRA{&chain, &res, nullptr};
            }
        }
    }
    ASSERT_NE(target_cra.residue, nullptr) << "No suitable residue found";

    // Case 1: applicable residue with apply_swap=true
    {
        int result = coot::rotamer_helper::optimize_rotamer_by_atom_names_gemmi(target_cra, true);
        EXPECT_EQ(result, 0);
    }

    // Case 2: same residue with apply_swap=false
    {
        int result = coot::rotamer_helper::optimize_rotamer_by_atom_names_gemmi(target_cra, false);
        EXPECT_EQ(result, 0);
    }

    // Case 3: non-applicable residue (ALA)
    {
        gemmi::CRA ala_cra{nullptr, nullptr, nullptr};
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.name == "ALA") {
                    ala_cra = gemmi::CRA{&chain, &res, nullptr};
                }
            }
        }
        if (ala_cra.residue) {
            int result = coot::rotamer_helper::optimize_rotamer_by_atom_names_gemmi(ala_cra, true);
            EXPECT_EQ(result, 0);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}