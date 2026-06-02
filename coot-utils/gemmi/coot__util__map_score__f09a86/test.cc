#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include "api/molecules-container.hh"
#include "api/coot-molecule.hh"
#include "coot-utils/coot-map-utils.hh"
#include "function.hh"

TEST(OracleTest, map_score) {
    // Setup: construct molecules container and load data
    molecules_container_t mc;
    mc.geometry_init_standard();

    int imol = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_GE(imol, 0);

    int imol_map = mc.read_mtz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
                               "FWT", "PHWT", "", false, false);
    EXPECT_GE(imol_map, 0);

    // Get the gemmi structure from the PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model &model = st.models[0];
    gemmi::Chain &chain = model.chains[0];

    // Get xmap from the map molecule
    const clipper::Xmap<float> &xmap = mc[imol_map].xmap;

    // Case 1: first residue (MET)
    {
        gemmi::Residue &res = chain.residues[0];

        EXPECT_EQ(res.atoms.size(), 8u);

        EXPECT_EQ(chain.name, "A");
        EXPECT_EQ(res.seqid.num.value, 1);
        EXPECT_EQ(res.name, "MET");

        std::vector<const gemmi::Atom *> atom_vec;
        for (gemmi::Atom &a : res.atoms) {
            atom_vec.push_back(&a);
        }

        float score = coot::util::map_score_gemmi(atom_vec, xmap);
        EXPECT_NEAR(score, 2.59989f, 1e-4);
    }

    // Case 2: empty atom vector (edge case)
    {
        std::vector<const gemmi::Atom *> empty_vec;
        EXPECT_EQ(empty_vec.size(), 0u);

        float empty_score = coot::util::map_score_gemmi(empty_vec, xmap);
        EXPECT_NEAR(empty_score, 0.0f, 1e-4);
    }

    // Case 3: different residue with more atoms (GLU)
    {
        gemmi::Residue &res2 = chain.residues[1];

        EXPECT_EQ(res2.atoms.size(), 9u);

        EXPECT_EQ(chain.name, "A");
        EXPECT_EQ(res2.seqid.num.value, 2);
        EXPECT_EQ(res2.name, "GLU");

        std::vector<const gemmi::Atom *> atom_vec2;
        for (gemmi::Atom &a : res2.atoms) {
            atom_vec2.push_back(&a);
        }

        float score2 = coot::util::map_score_gemmi(atom_vec2, xmap);
        EXPECT_NEAR(score2, 3.69131f, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}