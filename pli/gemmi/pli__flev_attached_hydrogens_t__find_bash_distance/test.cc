#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, find_bash_distance) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    gemmi::Model &model = st.models[0];

    std::vector<gemmi::CRA> close_atoms;
    if (!model.chains.empty()) {
        gemmi::Chain &chain = model.chains[0];
        if (!chain.residues.empty()) {
            gemmi::Residue &residue = chain.residues[0];
            int count = 0;
            for (gemmi::Atom &atom : residue.atoms) {
                close_atoms.push_back(gemmi::CRA{&chain, &residue, &atom});
                count++;
                if (count >= 5) break;
            }
        }
    }

    // Case 1: Place ligand_pos very close to a real atom so cannonball hits
    {
        double ref_x = close_atoms[0].atom->pos.x;
        double ref_y = close_atoms[0].atom->pos.y;
        double ref_z = close_atoms[0].atom->pos.z;

        clipper::Coord_orth hydrogen_pos(ref_x, ref_y, ref_z + 1.0);
        clipper::Coord_orth ligand_pos(ref_x, ref_y, ref_z + 2.0);

        coot::bash_distance_t result = pli::find_bash_distance_gemmi(ligand_pos, hydrogen_pos, close_atoms);

        EXPECT_NEAR(result.dist, 0.0, 1e-4);
        EXPECT_TRUE(result.limited);
    }

    // Case 2: Empty close_atoms — should return unlimited
    {
        std::vector<gemmi::CRA> empty_atoms;
        clipper::Coord_orth ligand_pos(10.0, 10.0, 10.0);
        clipper::Coord_orth hydrogen_pos(10.0, 10.0, 11.0);

        coot::bash_distance_t result = pli::find_bash_distance_gemmi(ligand_pos, hydrogen_pos, empty_atoms);

        EXPECT_NEAR(result.dist, -1.0, 1e-4);
        EXPECT_FALSE(result.limited);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}