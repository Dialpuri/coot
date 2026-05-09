#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include "function.hh"

// Helper function to compute center of mass (matches oracle logic)
static coot::Cartesian compute_center(const gemmi::Structure& st) {
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    int n_atoms = 0;
    
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            for (const auto& residue : chain.residues) {
                for (const auto& atom : residue.atoms) {
                    sum_x += atom.pos.x;
                    sum_y += atom.pos.y;
                    sum_z += atom.pos.z;
                    n_atoms++;
                }
            }
        }
    }
    
    if (n_atoms == 0) {
        return coot::Cartesian(0, 0, 0);
    }
    
    return coot::Cartesian(sum_x / n_atoms, sum_y / n_atoms, sum_z / n_atoms);
}

TEST(OracleTest, move_molecule_to_new_centre_basic_translation) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    coot::Cartesian current_center = compute_center(st);
    
    // Move by (10, 10, 10)
    coot::Cartesian new_centre(current_center.x() + 10.0,
                               current_center.y() + 10.0,
                               current_center.z() + 10.0);
    
    int status = coot::move_molecule_to_new_centre_gemmi(st, new_centre);
    EXPECT_EQ(status, 0);
    
    // Re-compute center after translation
    coot::Cartesian new_center_after = compute_center(st);
    
    // Verify new center matches target
    EXPECT_NEAR(new_center_after.x(), new_centre.x(), 1e-3);
    EXPECT_NEAR(new_center_after.y(), new_centre.y(), 1e-3);
    EXPECT_NEAR(new_center_after.z(), new_centre.z(), 1e-3);
    
    // Verify translation delta is (10, 10, 10)
    double delta_x = new_center_after.x() - current_center.x();
    double delta_y = new_center_after.y() - current_center.y();
    double delta_z = new_center_after.z() - current_center.z();
    
    EXPECT_NEAR(delta_x, 10.0, 1e-3);
    EXPECT_NEAR(delta_y, 10.0, 1e-3);
    EXPECT_NEAR(delta_z, 10.0, 1e-3);
}

TEST(OracleTest, move_molecule_to_new_centre_zero_translation) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    coot::Cartesian current_center = compute_center(st);
    
    // Move to current center (no translation)
    coot::Cartesian new_centre(current_center.x(), current_center.y(), current_center.z());
    
    int status = coot::move_molecule_to_new_centre_gemmi(st, new_centre);
    EXPECT_EQ(status, 0);
    
    // Re-compute center after "translation"
    coot::Cartesian new_center_after = compute_center(st);
    
    // Verify center hasn't changed
    EXPECT_NEAR(new_center_after.x(), current_center.x(), 1e-3);
    EXPECT_NEAR(new_center_after.y(), current_center.y(), 1e-3);
    EXPECT_NEAR(new_center_after.z(), current_center.z(), 1e-3);
    
    // Verify delta is (0, 0, 0)
    double delta_x = new_center_after.x() - current_center.x();
    double delta_y = new_center_after.y() - current_center.y();
    double delta_z = new_center_after.z() - current_center.z();
    
    EXPECT_NEAR(delta_x, 0.0, 1e-3);
    EXPECT_NEAR(delta_y, 0.0, 1e-3);
    EXPECT_NEAR(delta_z, 0.0, 1e-3);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}