#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

// Trim whitespace from both ends
static std::string trim(const std::string &s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

TEST(OracleTest, get_neighb_normal_gemmi) {
    // case: C=O bond normal calculation
    {
        // Read structure with gemmi
        gemmi::Structure st = gemmi::read_pdb_file(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        // Get the first model, first chain, first residue
        gemmi::Model* model = &st.models[0];
        gemmi::Chain* chain = &model->chains[0];
        gemmi::Residue* res = &chain->residues[0];

        // Copy atoms into a vector for indexed access
        std::vector<gemmi::Atom> atoms;
        for (const auto& a : res->atoms) {
            atoms.push_back(a);
        }
        int count = static_cast<int>(atoms.size());

        // Find backbone C and O atoms (trimming padded names)
        int iat_c = -1, iat_o = -1;
        for (int i = 0; i < count; i++) {
            std::string name = trim(atoms[i].name);
            if (name == "C") iat_c = i;
            if (name == "O") iat_o = i;
        }

        ASSERT_GE(iat_c, 0);
        ASSERT_GE(iat_o, 0);

        int imol = 0; // molecule index (not used in gemmi version)
        std::string res_name = res->name;

        coot::protein_geometry geom;
        geom.init_standard();

        // Case 1: C=O bond with also_2nd_order = false
        {
            clipper::Coord_orth normal = get_neighb_normal_gemmi(
                geom, imol, iat_c, iat_o, atoms, res_name, false);
            EXPECT_NEAR(normal.x(), -0.240966, 1e-4);
            EXPECT_NEAR(normal.y(), 0.942706, 1e-4);
            EXPECT_NEAR(normal.z(), -0.230739, 1e-4);
        }

        // Case 2: C=O bond with also_2nd_order = true
        {
            clipper::Coord_orth normal2 = get_neighb_normal_gemmi(
                geom, imol, iat_c, iat_o, atoms, res_name, true);
            EXPECT_NEAR(normal2.x(), -0.62607, 1e-4);
            EXPECT_NEAR(normal2.y(), 0.472419, 1e-4);
            EXPECT_NEAR(normal2.z(), 0.620368, 1e-4);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
