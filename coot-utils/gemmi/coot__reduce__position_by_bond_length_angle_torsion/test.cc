#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <clipper/core/clipper_types.h>
#include <clipper/core/coords.h>

TEST(OracleTest, position_by_bond_length_angle_torsion) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model *model = &st.models[0];
    ASSERT_NE(model, nullptr);

    const gemmi::Chain *chain = nullptr;
    for (const gemmi::Chain &c : model->chains) {
        if (c.name == "A") {
            chain = &c;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    const gemmi::Residue *res0 = nullptr;
    size_t idx = 0;
    for (const gemmi::Residue &r : chain->residues) {
        if (idx == 0) {
            res0 = &r;
            break;
        }
        ++idx;
    }
    ASSERT_NE(res0, nullptr);

    // gemmi atom.name is trimmed (e.g. "CA" not " CA ")
    const gemmi::Atom *ca = nullptr;
    const gemmi::Atom *cb = nullptr;
    for (const gemmi::Atom &a : res0->atoms) {
        if (a.name == "CA") ca = &a;
        if (a.name == "CB") cb = &a;
    }
    ASSERT_NE(ca, nullptr);
    ASSERT_NE(cb, nullptr);

    const gemmi::Atom *cg = nullptr;
    idx = 0;
    for (const gemmi::Residue &r : chain->residues) {
        if (idx == 1) {
            for (const gemmi::Atom &a : r.atoms) {
                if (a.name == "CG") {
                    cg = &a;
                    break;
                }
            }
            break;
        }
        ++idx;
    }
    ASSERT_NE(cg, nullptr);

    // Case 1: Typical call with CA, CB, CG atoms
    {
        double bl = 1.5;
        double angle_rad = clipper::Util::d2rad(109.5);
        double torsion_rad = clipper::Util::d2rad(60.0);

        clipper::Coord_orth pos = coot::reduce::position_by_bond_length_angle_torsion_gemmi(
            ca, cb, cg, bl, angle_rad, torsion_rad);

        EXPECT_NEAR(pos.x(), 13.3816, 1e-4);
        EXPECT_NEAR(pos.y(), 3.23037, 1e-4);
        EXPECT_NEAR(pos.z(), 80.7018, 1e-4);
    }

    // Case 2: Same atoms, different torsion angle (0 radians)
    {
        double bl = 1.09;
        double angle_rad = clipper::Util::d2rad(109.5);
        double torsion_rad = 0.0;

        clipper::Coord_orth pos = coot::reduce::position_by_bond_length_angle_torsion_gemmi(
            ca, cb, cg, bl, angle_rad, torsion_rad);

        EXPECT_NEAR(pos.x(), 12.5955, 1e-4);
        EXPECT_NEAR(pos.y(), 2.35747, 1e-4);
        EXPECT_NEAR(pos.z(), 80.2108, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}