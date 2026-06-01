#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, set_mol_cell) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GT(static_cast<int>(st.models.size()), 0);
    ASSERT_GT(static_cast<int>(st.models[0].chains.size()), 0);

    // Read original cell parameters before setting
    double orig_a = st.cell.a;
    double orig_b = st.cell.b;
    double orig_c = st.cell.c;
    double orig_alpha = st.cell.alpha;
    double orig_beta = st.cell.beta;
    double orig_gamma = st.cell.gamma;

    EXPECT_NEAR(orig_a, 53.512, 1e-4);
    EXPECT_NEAR(orig_b, 71.311, 1e-4);
    EXPECT_NEAR(orig_c, 72.101, 1e-4);
    EXPECT_NEAR(orig_alpha, 90, 1e-4);
    EXPECT_NEAR(orig_beta, 90, 1e-4);
    EXPECT_NEAR(orig_gamma, 90, 1e-4);

    // Case 1: set a reasonable orthorhombic cell (50/60/70 Angstroms, 90/90/90 degrees)
    {
        clipper::Cell_descr cd(50.0, 60.0, 70.0, 90.0, 90.0, 90.0);
        clipper::Cell cell_local(cd);

        EXPECT_NEAR(cell_local.a(), 50, 1e-4);
        EXPECT_NEAR(cell_local.b(), 60, 1e-4);
        EXPECT_NEAR(cell_local.c(), 70, 1e-4);
        EXPECT_NEAR(clipper::Util::rad2d(cell_local.alpha()), 90, 1e-4);
        EXPECT_NEAR(clipper::Util::rad2d(cell_local.beta()), 90, 1e-4);
        EXPECT_NEAR(clipper::Util::rad2d(cell_local.gamma()), 90, 1e-4);

        bool result = coot::util::set_mol_cell_gemmi(st, cell_local);
        EXPECT_TRUE(result);

        double new_a = st.cell.a;
        double new_b = st.cell.b;
        double new_c = st.cell.c;
        double new_alpha = st.cell.alpha;
        double new_beta = st.cell.beta;
        double new_gamma = st.cell.gamma;
        EXPECT_NEAR(new_a, 50, 1e-4);
        EXPECT_NEAR(new_b, 60, 1e-4);
        EXPECT_NEAR(new_c, 70, 1e-4);
        EXPECT_NEAR(new_alpha, 90, 1e-4);
        EXPECT_NEAR(new_beta, 90, 1e-4);
        EXPECT_NEAR(new_gamma, 90, 1e-4);
    }

    // Case 2: set a non-orthorhombic cell (monoclinic: alpha=90, beta=110, gamma=90)
    {
        clipper::Cell_descr cd2(40.0, 55.0, 65.0, 90.0, 110.0, 90.0);
        clipper::Cell cell_local2(cd2);

        EXPECT_NEAR(cell_local2.a(), 40, 1e-4);
        EXPECT_NEAR(cell_local2.b(), 55, 1e-4);
        EXPECT_NEAR(cell_local2.c(), 65, 1e-4);
        EXPECT_NEAR(clipper::Util::rad2d(cell_local2.alpha()), 90, 1e-4);
        EXPECT_NEAR(clipper::Util::rad2d(cell_local2.beta()), 110, 1e-4);
        EXPECT_NEAR(clipper::Util::rad2d(cell_local2.gamma()), 90, 1e-4);

        bool result2 = coot::util::set_mol_cell_gemmi(st, cell_local2);
        EXPECT_TRUE(result2);

        double new_a2 = st.cell.a;
        double new_b2 = st.cell.b;
        double new_c2 = st.cell.c;
        double new_alpha2 = st.cell.alpha;
        double new_beta2 = st.cell.beta;
        double new_gamma2 = st.cell.gamma;
        EXPECT_NEAR(new_a2, 40, 1e-4);
        EXPECT_NEAR(new_b2, 55, 1e-4);
        EXPECT_NEAR(new_c2, 65, 1e-4);
        EXPECT_NEAR(new_alpha2, 90, 1e-4);
        EXPECT_NEAR(new_beta2, 110, 1e-4);
        EXPECT_NEAR(new_gamma2, 90, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}