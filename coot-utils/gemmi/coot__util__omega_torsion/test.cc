#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: read a PDB file into a gemmi::Structure
static gemmi::Structure read_pdb(const std::string& file_name) {
    return gemmi::read_pdb_file(file_name);
}

// Test omega_torsion function
TEST(OracleTest, OmegaTorsion) {
    // case: normal residues
    {
        std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
        gemmi::Structure st = read_pdb(pdb_path);

        ASSERT_FALSE(st.models.empty());
        gemmi::Model& model = st.models[0];

        ASSERT_FALSE(model.chains.empty());
        gemmi::Chain& chain = model.chains[0];

        ASSERT_GE(chain.residues.size(), 2);

        const gemmi::Residue& res1 = chain.residues[0];
        const gemmi::Residue& res2 = chain.residues[1];

        std::string altconf("");
        std::pair<bool, double> result = coot::util::omega_torsion_gemmi(&res2, &res1, altconf);

        EXPECT_TRUE(result.first);
        EXPECT_NEAR(result.second, 1.45094, 1e-4);
        EXPECT_NEAR(result.second * 180.0 / M_PI, 83.1329, 1e-4);
    }

    // case: null C_residue
    {
        std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
        gemmi::Structure st = read_pdb(pdb_path);

        ASSERT_FALSE(st.models.empty());
        gemmi::Model& model = st.models[0];

        ASSERT_FALSE(model.chains.empty());
        gemmi::Chain& chain = model.chains[0];

        ASSERT_GE(chain.residues.size(), 10);
        const gemmi::Residue& res1 = chain.residues[9]; // //A/10

        std::string altconf("");
        std::pair<bool, double> r = coot::util::omega_torsion_gemmi(nullptr, &res1, altconf);

        EXPECT_FALSE(r.first);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}