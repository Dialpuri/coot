#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test for get_residue_centre_gemmi
TEST(OracleTest, GetResidueCentre) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case 1: valid residue with atoms (ILE 10)
    {
        bool found = false;
        for (gemmi::Model& model : st.models)
        for (gemmi::Chain& chain : model.chains)
        for (gemmi::Residue& res : chain.residues) {
            if (res.name == "ILE" && res.seqid.num.value == 10) {
                std::pair<bool, clipper::Coord_orth> result = coot::util::get_residue_centre_gemmi(res);
                EXPECT_TRUE(result.first) << "get_residue_centre_gemmi status should be true";
                if (result.first) {
                    clipper::Coord_orth c = result.second;
                    EXPECT_NEAR(c.x(), 27.2126, 1e-4);
                    EXPECT_NEAR(c.y(), 11.192, 1e-4);
                    EXPECT_NEAR(c.z(), 67.5274, 1e-4);
                }
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found) << "ILE 10 not found";
    }

    // case 2: valid residue with multiple atoms (LYS 20)
    {
        bool found = false;
        for (gemmi::Model& model : st.models)
        for (gemmi::Chain& chain : model.chains)
        for (gemmi::Residue& res : chain.residues) {
            if (res.name == "LYS" && res.seqid.num.value == 20) {
                std::pair<bool, clipper::Coord_orth> result = coot::util::get_residue_centre_gemmi(res);
                EXPECT_TRUE(result.first) << "get_residue_centre_gemmi status should be true";
                if (result.first) {
                    clipper::Coord_orth c = result.second;
                    EXPECT_NEAR(c.x(), 25.6588, 1e-4);
                    EXPECT_NEAR(c.y(), 8.0696, 1e-4);
                    EXPECT_NEAR(c.z(), 72.2372, 1e-4);
                }
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found) << "LYS 20 not found";
    }

    // case 3: residue with only 1 atom (HOH 1 in chain B)
    {
        bool found_single_atom_res = false;
        for (gemmi::Model& model : st.models)
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.atoms.size() == 1) {
                    std::pair<bool, clipper::Coord_orth> result = coot::util::get_residue_centre_gemmi(res);
                    EXPECT_TRUE(result.first) << "get_residue_centre_gemmi status should be true for 1-atom residue";
                    if (result.first) {
                        clipper::Coord_orth c = result.second;
                        EXPECT_NEAR(c.x(), 30.744, 1e-4);
                        EXPECT_NEAR(c.y(), 8.475, 1e-4);
                        EXPECT_NEAR(c.z(), 36.748, 1e-4);
                    }
                    found_single_atom_res = true;
                    break;
                }
            }
            if (found_single_atom_res) break;
        }
        EXPECT_TRUE(found_single_atom_res) << "No single-atom residue found in the structure";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}