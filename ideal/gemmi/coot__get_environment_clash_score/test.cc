#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetEnvironmentClashScore) {
    // Load PDB file using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get the first model
    const gemmi::Model& model = st.models[0];
    
    // Collect all atoms from the first model
    std::vector<const gemmi::Atom*> all_atoms;
    for (const gemmi::Model& m : st.models) {
        for (const gemmi::Chain& chain : m.chains) {
            for (const gemmi::Residue& res : chain.residues) {
                for (const gemmi::Atom& atom : res.atoms) {
                    all_atoms.push_back(&atom);
                }
            }
        }
    }
    
    // Test case 1: empty avoid list - should return 0
    {
        std::vector<std::pair<bool, clipper::Coord_orth>> avoid_these_atoms;
        double score = coot::get_environment_clash_score_gemmi(model, all_atoms, avoid_these_atoms);
        EXPECT_NEAR(score, 0, 1e-4);
    }
    
    // Test case 2: add an atom very close to the first atom in the selection
    // Get the first atom's coordinates
    if (!all_atoms.empty()) {
        clipper::Coord_orth close_atom = coot::co_gemmi(all_atoms[0]);
        // Add an atom very close (within 2.5 angstroms) - this should trigger a clash
        std::vector<std::pair<bool, clipper::Coord_orth>> avoid_these_atoms;
        // Make it 2.0 angstroms away - within the 2.5 limit for non-waters
        clipper::Coord_orth nearby(close_atom.x() + 2.0, close_atom.y(), close_atom.z());
        avoid_these_atoms.push_back(std::make_pair(false, nearby)); // false = not a water
        
        double score = coot::get_environment_clash_score_gemmi(model, all_atoms, avoid_these_atoms);
        EXPECT_NEAR(score, 8.77909, 1e-4);
    }
    
    // Test case 3: add an atom far away - should return 0
    {
        std::vector<std::pair<bool, clipper::Coord_orth>> avoid_these_atoms;
        clipper::Coord_orth far_atom(100.0, 100.0, 100.0);
        avoid_these_atoms.push_back(std::make_pair(false, far_atom));
        
        double score = coot::get_environment_clash_score_gemmi(model, all_atoms, avoid_these_atoms);
        EXPECT_NEAR(score, 0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}