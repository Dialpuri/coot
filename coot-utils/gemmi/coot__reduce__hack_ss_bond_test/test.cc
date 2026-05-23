#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, hack_ss_bond_test) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "load failed";

    const gemmi::Model& model = st.models[0];

    // case: CYS 118 with SG atom
    {
        const gemmi::Chain* chain_a = nullptr;
        for (const auto& c : model.chains) {
            if (c.name == "A") {
                chain_a = &c;
                break;
            }
        }
        ASSERT_NE(chain_a, nullptr) << "Chain A not found";
        
        const gemmi::Residue* cys_res = nullptr;
        for (const auto& r : chain_a->residues) {
            if (r.seqid.num.value == 118 && r.name == "CYS") {
                cys_res = &r;
                break;
            }
        }
        ASSERT_NE(cys_res, nullptr) << "CYS 118 not found";
        
        bool found_sg = false;
        for (const auto& at : cys_res->atoms) {
            if (at.name == "SG") {
                found_sg = true;
                break;
            }
        }
        ASSERT_TRUE(found_sg) << "CYS 118 should have SG atom";

        bool result = coot::reduce::hack_ss_bond_test_gemmi(cys_res, model);
        EXPECT_FALSE(result) << "CYS 118 hack_ss_bond_test should return false";
    }

    // case: CYS 177 with SG atom
    {
        const gemmi::Chain* chain_a = nullptr;
        for (const auto& c : model.chains) {
            if (c.name == "A") {
                chain_a = &c;
                break;
            }
        }
        ASSERT_NE(chain_a, nullptr) << "Chain A not found";
        
        const gemmi::Residue* cys_res2 = nullptr;
        for (const auto& r : chain_a->residues) {
            if (r.seqid.num.value == 177 && r.name == "CYS") {
                cys_res2 = &r;
                break;
            }
        }
        ASSERT_NE(cys_res2, nullptr) << "CYS 177 not found";
        
        bool found_sg = false;
        for (const auto& at : cys_res2->atoms) {
            if (at.name == "SG") {
                found_sg = true;
                break;
            }
        }
        ASSERT_TRUE(found_sg) << "CYS 177 should have SG atom";

        bool result = coot::reduce::hack_ss_bond_test_gemmi(cys_res2, model);
        EXPECT_FALSE(result) << "CYS 177 hack_ss_bond_test should return false";
    }

    // case: Non-CYS residue (ILE 10)
    {
        const gemmi::Chain* chain_a = nullptr;
        for (const auto& c : model.chains) {
            if (c.name == "A") {
                chain_a = &c;
                break;
            }
        }
        ASSERT_NE(chain_a, nullptr) << "Chain A not found";
        
        const gemmi::Residue* ile_res = nullptr;
        for (const auto& r : chain_a->residues) {
            if (r.seqid.num.value == 10 && r.name == "ILE") {
                ile_res = &r;
                break;
            }
        }
        ASSERT_NE(ile_res, nullptr) << "ILE 10 not found";

        bool result = coot::reduce::hack_ss_bond_test_gemmi(ile_res, model);
        EXPECT_FALSE(result) << "Non-CYS residue should return false";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}