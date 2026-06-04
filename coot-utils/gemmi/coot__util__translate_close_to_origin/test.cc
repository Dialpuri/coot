#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include "function.hh"

// Pad atom name to 4 chars like MMDB's GetAtomName (leading space for short names)
static std::string pad4(const std::string& n) {
    std::string s = n.size() <= 2 ? (" " + n) : n;
    s.resize(4, ' ');
    return s;
}

TEST(OracleTestGemmi, translate_close_to_origin) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: Record centroid before translation
    {
        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];
        gemmi::Residue& res = chain.residues[0];
        gemmi::Atom& atom = res.atoms[0];

        double x_before = atom.pos.x;
        double y_before = atom.pos.y;
        double z_before = atom.pos.z;

        std::string atom_label = pad4(atom.name);
        EXPECT_EQ(atom_label, " N  ");
        EXPECT_NEAR(x_before, 12.334, 1e-4);
        EXPECT_NEAR(y_before, 2.772, 1e-4);
        EXPECT_NEAR(z_before, 74.951, 1e-4);
        EXPECT_EQ(gemmi::count_atom_sites(st), 2107);

        // Call the function
        coot::util::translate_close_to_origin_gemmi(st);

        double x_after = atom.pos.x;
        double y_after = atom.pos.y;
        double z_after = atom.pos.z;

        EXPECT_NEAR(x_after, -41.178, 1e-4);
        EXPECT_NEAR(y_after, 2.772, 1e-4);
        EXPECT_NEAR(z_after, 2.85, 1e-4);
        EXPECT_NEAR(x_after - x_before, -53.512, 1e-4);
        EXPECT_NEAR(y_after - y_before, 0.0, 1e-4);
        EXPECT_NEAR(z_after - z_before, -72.101, 1e-4);
    }

    // Case 2: Try again on already-shifted molecule — shift should be zero
    {
        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];
        gemmi::Residue& res = chain.residues[0];
        gemmi::Atom& atom = res.atoms[0];

        double x_before2 = atom.pos.x;
        double y_before2 = atom.pos.y;
        double z_before2 = atom.pos.z;

        EXPECT_NEAR(x_before2, -41.178, 1e-4);
        EXPECT_NEAR(y_before2, 2.772, 1e-4);
        EXPECT_NEAR(z_before2, 2.85, 1e-4);

        coot::util::translate_close_to_origin_gemmi(st);

        double x_after2 = atom.pos.x;
        double y_after2 = atom.pos.y;
        double z_after2 = atom.pos.z;

        EXPECT_NEAR(x_after2, -41.178, 1e-4);
        EXPECT_NEAR(y_after2, 2.772, 1e-4);
        EXPECT_NEAR(z_after2, 2.85, 1e-4);
        EXPECT_NEAR(x_after2 - x_before2, 0.0, 1e-4);
        EXPECT_NEAR(y_after2 - y_before2, 0.0, 1e-4);
        EXPECT_NEAR(z_after2 - z_before2, 0.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}