#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include "function.hh"

TEST(OracleTest, ChiTorsionCalculation) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case 1: ASP residue (has chi angles)
    {
        const gemmi::Model& model = st.models[0];
        const gemmi::Chain& chain = model.chains[0];
        
        // Find ASP residue 68
        const gemmi::Residue* asp_res = nullptr;
        for (const auto& res : chain.residues) {
            if (res.seqid.num.value == 68 && res.name == "ASP") {
                asp_res = &res;
                break;
            }
        }
        ASSERT_NE(asp_res, nullptr) << "ASP residue not found";
        
        // ASSERT the residue name
        EXPECT_STREQ(asp_res->name.c_str(), "ASP");
        
        // Get the 4 chi angle atoms for ASP (N, CA, CB, CG)
        std::vector<int> atom_indices;
        std::vector<std::string> atom_names = {"N", "CA", "CB", "CG"};
        for (const auto& name : atom_names) {
            int idx = -1;
            for (size_t i = 0; i < asp_res->atoms.size(); ++i) {
                if (asp_res->atoms[i].name == name) {
                    idx = static_cast<int>(i);
                    break;
                }
            }
            atom_indices.push_back(idx);
        }
        
        // Call the gemmi ported function
        double chi = coot::rotamer::chi_torsion_gemmi(atom_indices, asp_res->atoms);
        EXPECT_NEAR(chi, -64.7521, 1e-4);
    }
    
    // case 2: GLY residue (no side chain chi angles)
    {
        const gemmi::Model& model = st.models[0];
        const gemmi::Chain& chain = model.chains[0];
        
        // Find GLY residue 11
        const gemmi::Residue* gly_res = nullptr;
        for (const auto& res : chain.residues) {
            if (res.seqid.num.value == 11 && res.name == "GLY") {
                gly_res = &res;
                break;
            }
        }
        ASSERT_NE(gly_res, nullptr) << "GLY residue not found";
        
        // ASSERT the residue name
        EXPECT_STREQ(gly_res->name.c_str(), "GLY");
        
        // Try to get chi angle atoms for GLY (should be empty or not have 4 atoms)
        std::vector<int> atom_indices;
        std::vector<std::string> atom_names = {"N", "CA", "CB", "CG"};
        for (const auto& name : atom_names) {
            int idx = -1;
            for (size_t i = 0; i < gly_res->atoms.size(); ++i) {
                if (gly_res->atoms[i].name == name) {
                    idx = static_cast<int>(i);
                    break;
                }
            }
            // CB doesn't exist in GLY, so we stop
            if (idx < 0) break;
            atom_indices.push_back(idx);
        }
        
        // GLY should not have 4 chi angle atoms (CB doesn't exist)
        EXPECT_TRUE(!(atom_indices.size() == 4))
            << "GLY should not have chi angles with 4 atoms";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}