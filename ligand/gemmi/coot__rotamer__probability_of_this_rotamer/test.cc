#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(RotamerProbabilityTest, ValidResidueWithSideChain) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];
    
    // Case 1: Valid residue with side chain (ILE 10 in chain A)
    gemmi::Chain* chain_a = nullptr;
    gemmi::Residue* res_ile = nullptr;
    for (auto& c : model.chains) {
        if (c.name == "A") {
            chain_a = &c;
            for (auto& r : c.residues) {
                if (r.seqid.num.value == 10 && r.name == "ILE") {
                    res_ile = &r;
                    break;
                }
            }
            break;
        }
    }
    
    ASSERT_NE(chain_a, nullptr);
    ASSERT_NE(res_ile, nullptr);
    
    gemmi::CRA cra{chain_a, res_ile, nullptr};
    coot::rotamer_probability_info_t info = 
        coot::rotamer_helper::probability_of_this_rotamer_gemmi(cra);
    
    EXPECT_EQ(info.state, 1);
    EXPECT_NEAR(info.probability, 13, 1e-4);
    EXPECT_EQ(info.rotamer_name, "pt");
}

TEST(RotamerProbabilityTest, GLYResidueNoSideChain) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];
    
    // Case 2: GLY residue (no side chain, should return state=-2)
    gemmi::Chain* chain_a = nullptr;
    gemmi::Residue* res_gly = nullptr;
    for (auto& c : model.chains) {
        if (c.name == "A") {
            chain_a = &c;
            for (auto& r : c.residues) {
                if (r.seqid.num.value == 11 && r.name == "GLY") {
                    res_gly = &r;
                    break;
                }
            }
            break;
        }
    }
    
    ASSERT_NE(chain_a, nullptr);
    ASSERT_NE(res_gly, nullptr);
    
    gemmi::CRA cra{chain_a, res_gly, nullptr};
    coot::rotamer_probability_info_t info = 
        coot::rotamer_helper::probability_of_this_rotamer_gemmi(cra);
    
    EXPECT_EQ(info.state, -2);
    EXPECT_NEAR(info.probability, 0, 1e-4);
    EXPECT_EQ(info.rotamer_name, "none");
}

TEST(RotamerProbabilityTest, ALAResidueNoSideChain) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];
    
    // Case 3: ALA residue (no side chain, should return state=-2)
    gemmi::Chain* chain_a = nullptr;
    gemmi::Residue* res_ala = nullptr;
    for (auto& c : model.chains) {
        if (c.name == "A") {
            chain_a = &c;
            for (auto& r : c.residues) {
                if (r.seqid.num.value == 21 && r.name == "ALA") {
                    res_ala = &r;
                    break;
                }
            }
            break;
        }
    }
    
    ASSERT_NE(chain_a, nullptr);
    ASSERT_NE(res_ala, nullptr);
    
    gemmi::CRA cra{chain_a, res_ala, nullptr};
    coot::rotamer_probability_info_t info = 
        coot::rotamer_helper::probability_of_this_rotamer_gemmi(cra);
    
    EXPECT_EQ(info.state, -2);
    EXPECT_NEAR(info.probability, 0, 1e-4);
    EXPECT_EQ(info.rotamer_name, "none");
}

TEST(RotamerProbabilityTest, ValidResiduePHE) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];
    
    // Case 4: Valid residue with side chain (PHE 4 in chain A)
    gemmi::Chain* chain_a = nullptr;
    gemmi::Residue* res_phe = nullptr;
    for (auto& c : model.chains) {
        if (c.name == "A") {
            chain_a = &c;
            for (auto& r : c.residues) {
                if (r.seqid.num.value == 4 && r.name == "PHE") {
                    res_phe = &r;
                    break;
                }
            }
            break;
        }
    }
    
    ASSERT_NE(chain_a, nullptr);
    ASSERT_NE(res_phe, nullptr);
    
    gemmi::CRA cra{chain_a, res_phe, nullptr};
    coot::rotamer_probability_info_t info = 
        coot::rotamer_helper::probability_of_this_rotamer_gemmi(cra);
    
    EXPECT_EQ(info.state, 1);
    EXPECT_NEAR(info.probability, 44, 1e-4);
    EXPECT_EQ(info.rotamer_name, "m-85");
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}