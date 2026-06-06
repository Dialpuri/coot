#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "geometry/protein-geometry.hh"
#include "mini-mol/atom-quads.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_residue/gemmi/function.hh"

TEST(OracleTest, rotate_around_bond_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: valid residue with torsion atoms (residue 2 GLU, quad N-CA-CB-CG)
    {
        std::string residue_cid = "//A/2";
        std::string alt_conf = "";
        coot::atom_name_quad quad("N", "CA", "CB", "CG");
        double torsion_angle = 90.0;

        int result = coot::rotate_around_bond_gemmi(residue_cid, alt_conf, quad, torsion_angle, geom, model);

        EXPECT_EQ(result, 0);

        // Verify residue was found
        gemmi::Residue* res = coot::cid_to_residue_gemmi(residue_cid, model);
        EXPECT_NE(res, nullptr);
        if (res) {
            EXPECT_STREQ(res->name.c_str(), "GLU");
            EXPECT_EQ(res->seqid.num.value, 2);
        }
    }

    // Case 2: invalid residue CID (non-existent residue number)
    {
        std::string residue_cid = "//A/99999";
        std::string alt_conf = "";
        coot::atom_name_quad quad("N", "CA", "CB", "CG");
        double torsion_angle = 0.0;

        int result = coot::rotate_around_bond_gemmi(residue_cid, alt_conf, quad, torsion_angle, geom, model);

        EXPECT_EQ(result, 0);

        gemmi::Residue* res = coot::cid_to_residue_gemmi(residue_cid, model);
        EXPECT_EQ(res, nullptr);
    }

    // Case 3: valid residue but with blank atom name (should fail inside set_dihedral)
    {
        std::string residue_cid = "//A/10";
        std::string alt_conf = "";
        coot::atom_name_quad quad("", "CA", "CB", "CG");
        double torsion_angle = 0.0;

        int result = coot::rotate_around_bond_gemmi(residue_cid, alt_conf, quad, torsion_angle, geom, model);

        EXPECT_EQ(result, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}