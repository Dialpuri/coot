#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, IsMainChainP) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model& model = st.models[0];
    
    auto find_atom = [&model](const std::string& chain_id, int resnum, const std::string& inscode,
                               const std::string& atom_name, const gemmi::Residue*& out_res) -> const gemmi::Atom* {
        for (auto& chain : model.chains) {
            if (chain.name == chain_id) {
                for (auto& residue : chain.residues) {
                    if (residue.seqid.num.value == resnum) {
                        std::string res_inscode(1, residue.seqid.icode);
                        if (res_inscode == inscode) {
                            out_res = &residue;
                            for (auto& atom : residue.atoms) {
                                if (atom.name == atom_name) {
                                    return &atom;
                                }
                            }
                        }
                    }
                }
            }
        }
        out_res = nullptr;
        return nullptr;
    };

    {
        const gemmi::Residue* res = nullptr;
        const gemmi::Atom* at = find_atom("A", 1, " ", "N", res);
        ASSERT_NE(at, nullptr);
        ASSERT_NE(res, nullptr);
        bool result = coot::is_main_chain_p_gemmi(*at, *res);
        EXPECT_TRUE(result);
    }

    {
        const gemmi::Residue* res = nullptr;
        const gemmi::Atom* at = find_atom("A", 1, " ", "CA", res);
        ASSERT_NE(at, nullptr);
        ASSERT_NE(res, nullptr);
        bool result = coot::is_main_chain_p_gemmi(*at, *res);
        EXPECT_TRUE(result);
    }

    {
        const gemmi::Residue* res = nullptr;
        const gemmi::Atom* at = find_atom("A", 1, " ", "CB", res);
        ASSERT_NE(at, nullptr);
        ASSERT_NE(res, nullptr);
        bool result = coot::is_main_chain_p_gemmi(*at, *res);
        EXPECT_FALSE(result);
    }

    {
        const gemmi::Residue* res = nullptr;
        const gemmi::Atom* at = find_atom("A", 11, " ", "HA2", res);
        ASSERT_NE(at, nullptr);
        ASSERT_NE(res, nullptr);
        bool result = coot::is_main_chain_p_gemmi(*at, *res);
        EXPECT_TRUE(result);
    }

    {
        const gemmi::Residue* res = nullptr;
        const gemmi::Atom* at = find_atom("A", 11, " ", "HA3", res);
        ASSERT_NE(at, nullptr);
        ASSERT_NE(res, nullptr);
        bool result = coot::is_main_chain_p_gemmi(*at, *res);
        EXPECT_TRUE(result);
    }

    {
        const gemmi::Residue* res = nullptr;
        const gemmi::Atom* at = find_atom("A", 298, " ", "OXT", res);
        ASSERT_NE(at, nullptr);
        ASSERT_NE(res, nullptr);
        bool result = coot::is_main_chain_p_gemmi(*at, *res);
        EXPECT_TRUE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}