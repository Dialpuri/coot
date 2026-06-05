#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <cmath>

TEST(OracleTest, position_by_bond_length_angle_torsion) {
    // Load PDB via gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];

    gemmi::Chain* chain = model.find_chain("A");
    ASSERT_NE(chain, nullptr);

    // Case 1: Valid atoms - ALA residue 10 chain A: CA, CB, C
    {
        gemmi::Residue* res = nullptr;
        for (gemmi::Residue& r : chain->residues) {
            if (r.seqid.num.value == 10) {
                res = &r;
                break;
            }
        }
        ASSERT_NE(res, nullptr);

        const gemmi::Atom* at_1 = res->find_atom("CA", '*');  // CA
        const gemmi::Atom* at_2 = res->find_atom("CB", '*');  // CB
        const gemmi::Atom* at_3 = res->find_atom("C",  '*');  // C

        ASSERT_NE(at_1, nullptr);
        ASSERT_NE(at_2, nullptr);
        ASSERT_NE(at_3, nullptr);

        double bl = 1.09;
        double angle_rad = 109.5 * M_PI / 180.0;
        double torsion_rad = 0.0;

        clipper::Coord_orth result = coot::reduce::position_by_bond_length_angle_torsion_gemmi(at_1, at_2, at_3, bl, angle_rad, torsion_rad);

        EXPECT_NEAR(result.x(), 26.6691, 1e-4);
        EXPECT_NEAR(result.y(), 13.5967, 1e-4);
        EXPECT_NEAR(result.z(), 66.903, 1e-4);
    }

    // Case 2: Different torsion angle (180 degrees)
    {
        gemmi::Residue* res = nullptr;
        for (gemmi::Residue& r : chain->residues) {
            if (r.seqid.num.value == 10) {
                res = &r;
                break;
            }
        }
        ASSERT_NE(res, nullptr);

        const gemmi::Atom* at_1 = res->find_atom("CA", '*');  // CA
        const gemmi::Atom* at_2 = res->find_atom("CB", '*');  // CB
        const gemmi::Atom* at_3 = res->find_atom("C",  '*');  // C

        ASSERT_NE(at_1, nullptr);
        ASSERT_NE(at_2, nullptr);
        ASSERT_NE(at_3, nullptr);

        double bl = 1.09;
        double angle_rad = 109.5 * M_PI / 180.0;
        double torsion_rad = M_PI; // 180 degrees

        clipper::Coord_orth result = coot::reduce::position_by_bond_length_angle_torsion_gemmi(at_1, at_2, at_3, bl, angle_rad, torsion_rad);

        EXPECT_NEAR(result.x(), 27.01, 1e-4);
        EXPECT_NEAR(result.y(), 12.4384, 1e-4);
        EXPECT_NEAR(result.z(), 65.2402, 1e-4);
    }

    // Case 3: Different residue - residue 20 chain A: N, CA, C
    {
        gemmi::Residue* res2 = nullptr;
        for (gemmi::Residue& r : chain->residues) {
            if (r.seqid.num.value == 20) {
                res2 = &r;
                break;
            }
        }
        ASSERT_NE(res2, nullptr);

        const gemmi::Atom* at_n = res2->find_atom("N",  '*');
        const gemmi::Atom* at_ca = res2->find_atom("CA", '*');
        const gemmi::Atom* at_c = res2->find_atom("C",  '*');

        ASSERT_NE(at_n, nullptr);
        ASSERT_NE(at_ca, nullptr);
        ASSERT_NE(at_c, nullptr);

        double bl = 1.09;
        double angle_rad = 109.5 * M_PI / 180.0;
        double torsion_rad = M_PI / 2.0; // 90 degrees

        clipper::Coord_orth result = coot::reduce::position_by_bond_length_angle_torsion_gemmi(at_n, at_ca, at_c, bl, angle_rad, torsion_rad);

        EXPECT_NEAR(result.x(), 24.316, 1e-4);
        EXPECT_NEAR(result.y(), 6.6811, 1e-4);
        EXPECT_NEAR(result.z(), 72.7011, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}