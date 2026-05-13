#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ... TEST(...) blocks here ...

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, transform_mol_identity_preserves_coords) {
    // Oracle used absolute path directly
    std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
    
    gemmi::Structure st = gemmi::read_pdb_file(pdb_path);
    
    // Get initial coordinates of first two atoms
    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];
    gemmi::Residue& res = chain.residues[0];
    gemmi::Atom& at1 = res.atoms[0];
    gemmi::Atom& at2 = res.atoms[1];

    double x1_before = at1.pos.x, y1_before = at1.pos.y, z1_before = at1.pos.z;
    double x2_before = at2.pos.x, y2_before = at2.pos.y, z2_before = at2.pos.z;

    // Identity transformation
    clipper::RTop_orth rtop_identity(
        clipper::Mat33<double>(1,0,0,0,1,0,0,0,1),
        clipper::Vec3<double>(0,0,0)
    );
    
    coot::util::transform_mol_gemmi(&st, rtop_identity);

    gemmi::Atom& at1_after = res.atoms[0];
    gemmi::Atom& at2_after = res.atoms[1];

    double x1_after = at1_after.pos.x, y1_after = at1_after.pos.y, z1_after = at1_after.pos.z;
    double x2_after = at2_after.pos.x, y2_after = at2_after.pos.y, z2_after = at2_after.pos.z;

    // EXPECT_EQ outputs for debugging
    std::cout << "INPUT  before: atom1(" << x1_before << "," << y1_before << "," << z1_before << ")" << std::endl;
    std::cout << "INPUT  before: atom2(" << x2_before << "," << y2_before << "," << z2_before << ")" << std::endl;
    std::cout << "OUTPUT after: atom1(" << x1_after << "," << y1_after << "," << z1_after << ")" << std::endl;
    std::cout << "OUTPUT after: atom2(" << x2_after << "," << y2_after << "," << z2_after << ")" << std::endl;

    // Assertions using expected values from oracle output
    EXPECT_NEAR(x1_after, 12.334, 1e-4);
    EXPECT_NEAR(y1_after, 2.772, 1e-4);
    EXPECT_NEAR(z1_after, 74.951, 1e-4);
    EXPECT_NEAR(x2_after, 13.401, 1e-4);
    EXPECT_NEAR(y2_after, 3.78, 1e-4);
    EXPECT_NEAR(z2_after, 74.703, 1e-4);

    // Check if coordinates are unchanged after identity transform
    bool unchanged = (std::abs(x1_after - x1_before) < 1e-6 &&
                      std::abs(y1_after - y1_before) < 1e-6 &&
                      std::abs(z1_after - z1_before) < 1e-6 &&
                      std::abs(x2_after - x2_before) < 1e-6 &&
                      std::abs(y2_after - y2_before) < 1e-6 &&
                      std::abs(z2_after - z2_before) < 1e-6);
    EXPECT_TRUE(unchanged) << "identity_transform_preserves_coords should be true";
}