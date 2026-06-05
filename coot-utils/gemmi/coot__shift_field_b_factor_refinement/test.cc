#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/calculate.hpp>
#include "function.hh"

static double compute_b_sum(const gemmi::Structure& st) {
    double b_sum = 0;
    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& res : chain.residues) {
                for (const gemmi::Atom& atom : res.atoms) {
                    b_sum += atom.b_iso;
                }
            }
        }
    }
    return b_sum;
}

TEST(OracleTest, shift_field_b_factor_refinement) {
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif", gemmi::CoorFormat::Detect);
    ASSERT_FALSE(st.models.empty());

    double b_before = compute_b_sum(st);
    EXPECT_NEAR(b_before, 93060.9, 0.1);

    bool ok = coot::shift_field_b_factor_refinement_gemmi(st);
    // Refinement returns false because F_obs could not be loaded from MTZ
    EXPECT_FALSE(ok);

    double b_after = compute_b_sum(st);
    EXPECT_NEAR(b_after, 93060.9, 0.1);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}