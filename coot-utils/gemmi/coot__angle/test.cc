#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static const gemmi::Atom* get_atom_by_cid(const gemmi::Structure& st,
                                          const std::string& cid) {
    if (cid.size() < 5 || cid[0] != '/' || cid[1] != '/') {
        return nullptr;
    }
    
    size_t slash1 = cid.find('/', 2);
    if (slash1 == std::string::npos) return nullptr;
    
    size_t slash2 = cid.find('/', slash1 + 1);
    if (slash2 == std::string::npos) return nullptr;
    
    std::string chain_id = cid.substr(2, slash1 - 2);
    std::string res_str = cid.substr(slash1 + 1, slash2 - slash1 - 1);
    std::string atom_name = cid.substr(slash2 + 1);
    
    int res_num;
    try {
        res_num = std::stoi(res_str);
    } catch (...) {
        return nullptr;
    }
    
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name == chain_id) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == res_num) {
                        for (auto& atom : res.atoms) {
                            if (atom.name == atom_name) {
                                return &atom;
                            }
                        }
                    }
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, AngleFunction) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    {
        std::string cid_1 = "//A/10/CA";
        std::string cid_2 = "//A/10/C";
        std::string cid_3 = "//A/11/N";
        
        const gemmi::Atom *at_1 = get_atom_by_cid(st, cid_1);
        const gemmi::Atom *at_2 = get_atom_by_cid(st, cid_2);
        const gemmi::Atom *at_3 = get_atom_by_cid(st, cid_3);
        
        ASSERT_NE(at_1, nullptr);
        ASSERT_NE(at_2, nullptr);
        ASSERT_NE(at_3, nullptr);
        
        double angle_deg = coot::angle_gemmi(at_1, at_2, at_3);
        EXPECT_NEAR(angle_deg, 117.165, 0.001);
    }
    
    {
        std::string cid_1 = "//A/10/CA";
        std::string cid_2 = "//A/9999/N";
        std::string cid_3 = "//A/11/CA";
        
        const gemmi::Atom *at_1 = get_atom_by_cid(st, cid_1);
        const gemmi::Atom *at_2 = get_atom_by_cid(st, cid_2);
        const gemmi::Atom *at_3 = get_atom_by_cid(st, cid_3);
        
        EXPECT_EQ(at_2, nullptr);
        
        double angle_deg = coot::angle_gemmi(at_1, at_2, at_3);
        EXPECT_EQ(angle_deg, -1);
    }
    
    {
        std::string cid_1 = "//A/10/CA";
        std::string cid_2 = "//A/10/CA";
        std::string cid_3 = "//A/10/CA";
        
        const gemmi::Atom *at_1 = get_atom_by_cid(st, cid_1);
        const gemmi::Atom *at_2 = get_atom_by_cid(st, cid_2);
        const gemmi::Atom *at_3 = get_atom_by_cid(st, cid_3);
        
        ASSERT_NE(at_1, nullptr);
        ASSERT_NE(at_2, nullptr);
        ASSERT_NE(at_3, nullptr);
        
        double angle_deg = coot::angle_gemmi(at_1, at_2, at_3);
        EXPECT_NE(angle_deg, angle_deg);
    }
    
    {
        std::string cid_1 = "//A/10/N";
        std::string cid_2 = "//A/10/CA";
        std::string cid_3 = "//A/10/C";
        
        const gemmi::Atom *at_1 = get_atom_by_cid(st, cid_1);
        const gemmi::Atom *at_2 = get_atom_by_cid(st, cid_2);
        const gemmi::Atom *at_3 = get_atom_by_cid(st, cid_3);
        
        ASSERT_NE(at_1, nullptr);
        ASSERT_NE(at_2, nullptr);
        ASSERT_NE(at_3, nullptr);
        
        double angle_deg = coot::angle_gemmi(at_1, at_2, at_3);
        EXPECT_NEAR(angle_deg, 111.337, 0.001);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}