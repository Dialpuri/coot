#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "api/molecules-container.hh"

TEST(OracleTest, sum_density_for_atoms_in_residue) {
    molecules_container_t mc;

    // Load protein model
    int imol = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(imol, 0);

    // Load map from MTZ
    int imol_map = mc.read_mtz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz", "FWT", "PHWT", "", false, false);
    ASSERT_GE(imol_map, 0);

    clipper::Xmap<float> xmap = mc.get_xmap(imol_map);

    // Load PDB with gemmi for the model
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Case 1: Valid residue with real atom names (CA from residue 10, chain A)
    {
        std::string cid = "//A/10";
        std::vector<std::string> atom_names = {" CA "};
        double v = coot::sum_density_for_atoms_in_residue_gemmi(model, cid, atom_names, xmap);
        EXPECT_NEAR(v, 0.616902, 1e-4);
    }

    // Case 2: Valid residue with multiple atom names (backbone atoms from residue 50, chain A)
    {
        std::string cid = "//A/50";
        std::vector<std::string> atom_names = {" N  ", " CA ", " C  ", " O  "};
        double v = coot::sum_density_for_atoms_in_residue_gemmi(model, cid, atom_names, xmap);
        EXPECT_NEAR(v, 2.20078, 1e-4);
    }

    // Case 3: Invalid residue (residue number 9999 doesn't exist) — should return 0.0
    {
        std::string cid = "//A/9999";
        std::vector<std::string> atom_names = {" CA "};
        double v = coot::sum_density_for_atoms_in_residue_gemmi(model, cid, atom_names, xmap);
        EXPECT_NEAR(v, 0.0, 1e-4);
    }

    // Case 4: Valid residue but with atom name that doesn't exist
    {
        std::string cid = "//A/50";
        std::vector<std::string> atom_names = {" ZZZ "};
        double v = coot::sum_density_for_atoms_in_residue_gemmi(model, cid, atom_names, xmap);
        EXPECT_NEAR(v, 0.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}