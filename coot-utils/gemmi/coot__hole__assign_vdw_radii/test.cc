#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"
#include <clipper/clipper.h>

TEST(OracleTest, hole) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    int n_models = st.models.size();
    int n_atoms_total = gemmi::count_atom_sites(st);

    EXPECT_EQ(n_models, 1);
    EXPECT_EQ(n_atoms_total, 2107);

    // Initialise protein geometry
    coot::protein_geometry geom;
    geom.init_standard();

    // Create hole object using constructor
    clipper::Coord_orth from_pt(10.0, 10.0, 10.0);
    clipper::Coord_orth to_pt(50.0, 50.0, 50.0);

    EXPECT_NEAR(from_pt.x(), 10.0, 1e-4);
    EXPECT_NEAR(from_pt.y(), 10.0, 1e-4);
    EXPECT_NEAR(from_pt.z(), 10.0, 1e-4);
    EXPECT_NEAR(to_pt.x(), 50.0, 1e-4);
    EXPECT_NEAR(to_pt.y(), 50.0, 1e-4);
    EXPECT_NEAR(to_pt.z(), 50.0, 1e-4);

    // Assign VDW radii using gemmi API
    coot::hole::vdw_radii_result radii_result =
        coot::hole::assign_vdw_radii_gemmi(st.models[0], geom);

    // Read radii from first few atoms to verify assignment
    gemmi::Chain& first_chain = st.models[0].chains[0];
    gemmi::Residue& first_res = first_chain.residues[0];
    gemmi::Atom& first_atom = first_res.atoms[0];

    EXPECT_EQ(coot::hole::pad_atom_name(first_atom.name), " N  ");
    EXPECT_EQ(std::string(first_res.name), "MET");
    EXPECT_NEAR(radii_result.radii[0][0][0], 1.55, 1e-4);

    // Read radii from all atoms to compute statistics
    double radii_sum = radii_result.sum();
    int radii_count = radii_result.count();

    EXPECT_EQ(radii_count, 2107);
    EXPECT_EQ(n_atoms_total, 2107);
    EXPECT_NEAR(radii_sum / radii_count, 1.64234, 1e-4);

    // Case: empty molecule (no atoms, tests the !mol / no-atoms branch)
    {
        gemmi::Model empty_model = gemmi::Model();
        // Do not load any file — model is empty

        coot::hole::vdw_radii_result empty_result =
            coot::hole::assign_vdw_radii_gemmi(empty_model, geom);
        // No atoms to assign radii to — no crash
        EXPECT_EQ(empty_result.count(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}