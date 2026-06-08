#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ... TEST(...) blocks here ...

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, match_torsions_get_torsion) {
    // Load structure with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get chain A
    const gemmi::Chain *chain = nullptr;
    for (const auto &c : st.models[0].chains) {
        if (c.name == "A") {
            chain = &c;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    // Find residue GLU 2 (reference residue) and ASN 3 (moving residue)
    const gemmi::Residue *res_a = nullptr; // seqid 2
    const gemmi::Residue *res_b = nullptr; // seqid 3
    for (const auto &res : chain->residues) {
        if (res.seqid.num.value == 2 && !res.is_water()) {
            res_a = &res;
        }
        if (res.seqid.num.value == 3 && !res.is_water()) {
            res_b = &res;
        }
    }
    ASSERT_NE(res_a, nullptr);
    ASSERT_NE(res_b, nullptr);

    // Case 1: REFERENCE_TORSION with valid quad (N, CA, CB, CG from GLU res_a)
    {
        coot::atom_name_quad quad("N", "CA", "CB", "CG");
        int torsion_type = coot::match_torsions::REFERENCE_TORSION;
        std::pair<bool, double> result =
            coot::match_torsions::get_torsion_gemmi(torsion_type, *res_a, *res_b, quad);
        EXPECT_TRUE(result.first);
        EXPECT_NEAR(result.second, 1.50709, 1e-4);
    }

    // Case 2: MOVING_TORSION with valid quad (N, CA, CB, OD1 from ASN res_b)
    {
        coot::atom_name_quad quad("N", "CA", "CB", "OD1");
        int torsion_type = coot::match_torsions::MOVING_TORSION;
        std::pair<bool, double> result =
            coot::match_torsions::get_torsion_gemmi(torsion_type, *res_a, *res_b, quad);
        EXPECT_TRUE(result.first);
        EXPECT_NEAR(result.second, -0.7212, 1e-4);
    }

    // Case 3: invalid quad (atom name that doesn't exist)
    {
        coot::atom_name_quad quad("N", "CA", "CB", "FAKE");
        int torsion_type = coot::match_torsions::REFERENCE_TORSION;
        std::pair<bool, double> result =
            coot::match_torsions::get_torsion_gemmi(torsion_type, *res_a, *res_b, quad);
        EXPECT_FALSE(result.first);
        EXPECT_NEAR(result.second, 0.0, 1e-4);
    }

    // Case 4: invalid torsion_type (default branch) — returns {false, 0}
    {
        coot::atom_name_quad quad("N", "CA", "CB", "CG");
        int torsion_type = 999; // not REFERENCE or MOVING
        std::pair<bool, double> result =
            coot::match_torsions::get_torsion_gemmi(torsion_type, *res_a, *res_b, quad);
        EXPECT_FALSE(result.first);
        EXPECT_NEAR(result.second, 0.0, 1e-4);
    }
}
