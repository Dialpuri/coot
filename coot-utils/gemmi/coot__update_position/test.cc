#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, update_position) {
    // case: basic update_position with offset +10 in each direction
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Find the atom at chain A, residue 10, atom CA
    gemmi::Atom* atom = nullptr;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name == "A") {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 10) {
                        for (auto& a : res.atoms) {
                            if (a.name == "CA") {
                                atom = &a;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    
    ASSERT_NE(atom, nullptr) << "atom not found";
    
    // Record original coordinates
    float orig_x = atom->pos.x;
    float orig_y = atom->pos.y;
    float orig_z = atom->pos.z;
    
    // New position - offset by +10 Angstroms in each direction
    clipper::Coord_orth new_pos(orig_x + 10.0, orig_y + 10.0, orig_z + 10.0);
    
    // Call update_position_gemmi
    coot::update_position_gemmi(atom, new_pos);
    
    // Verify coordinates changed
    EXPECT_NEAR(atom->pos.x, 37.175, 1e-4);
    EXPECT_NEAR(atom->pos.y, 22.237, 1e-4);
    EXPECT_NEAR(atom->pos.z, 77.71, 1e-4);
    EXPECT_TRUE(atom->pos.x != orig_x || atom->pos.y != orig_y || atom->pos.z != orig_z);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}