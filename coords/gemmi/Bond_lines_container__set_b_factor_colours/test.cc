#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, set_b_factor_colours) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());
    gemmi::Chain& chain = model.chains[0];
    ASSERT_FALSE(chain.residues.empty());
    gemmi::Residue& res = chain.residues[0];
    ASSERT_FALSE(res.atoms.empty());
    gemmi::Atom& first_atom = res.atoms[0];

    Bond_lines_container bonds;

    // ── Case 1: standard protein PDB — normal call ──
    {
        bonds.b_factor_scale = 1.0f;
        int udd_handle_before = bonds.set_b_factor_colours_gemmi(st);

        int udd_handle = bonds.set_b_factor_colours_gemmi(st);

        float udd_val = -999;
        bool udd_read_ok = bonds.get_b_factor_fraction(chain, res, first_atom, udd_val);

        EXPECT_EQ(udd_handle, 16777217);
        EXPECT_EQ(udd_handle_before, 16777217);
        EXPECT_TRUE(udd_read_ok);
        EXPECT_NEAR(udd_val, 0.929286, 1e-4);
        EXPECT_EQ(first_atom.name, "N");
        EXPECT_NEAR(first_atom.b_iso, 65.05, 1e-4);
    }

    // ── Case 2: b_factor_scale = 0 → all fractions should be 0.0 ──
    {
        bonds.b_factor_scale = 0.0f;
        int udd_handle = bonds.set_b_factor_colours_gemmi(st);

        float udd_val = -999;
        bool udd_read_ok = bonds.get_b_factor_fraction(chain, res, first_atom, udd_val);

        EXPECT_EQ(udd_handle, 16777217);
        EXPECT_TRUE(udd_read_ok);
        EXPECT_NEAR(udd_val, 0.0, 1e-4);
    }

    // ── Case 3: large b_factor_scale → fractions should be clamped at 1.0 ──
    {
        bonds.b_factor_scale = 1000.0f;
        int udd_handle = bonds.set_b_factor_colours_gemmi(st);

        float udd_val = -999;
        bool udd_read_ok = bonds.get_b_factor_fraction(chain, res, first_atom, udd_val);

        EXPECT_EQ(udd_handle, 16777217);
        EXPECT_TRUE(udd_read_ok);
        EXPECT_NEAR(udd_val, 1.0, 1e-4);
        EXPECT_NEAR(first_atom.b_iso, 65.05, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}