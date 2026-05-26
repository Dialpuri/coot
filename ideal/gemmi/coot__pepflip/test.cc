#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static const gemmi::Residue* find_residue(const gemmi::Structure& st,
                                          const std::string& chain_id,
                                          int resno,
                                          const std::string& ins_code) {
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            if (chain.name == chain_id) {
                for (const auto& res : chain.residues) {
                    if (res.seqid.num.value == resno) {
                        std::string res_ins_code(1, res.seqid.icode);
                        if (res_ins_code == ins_code) {
                            return &res;
                        }
                    }
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, flip_peptide_function) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    {
        const gemmi::Residue* res = find_residue(st, "A", 10, " ");
        ASSERT_NE(res, nullptr) << "Case 1: residue 10 chain A not found";
        int result = coot::pepflip_gemmi(st, "A", 10, " ", "");
        EXPECT_EQ(result, 1) << "Case 1: flip_performed expected for CA atom on residue 10 chain A";
    }
    
    {
        const gemmi::Residue* res = find_residue(st, "A", 11, " ");
        ASSERT_NE(res, nullptr) << "Case 2: residue 11 chain A not found";
        int result = coot::pepflip_gemmi(st, "A", 11, " ", "");
        EXPECT_EQ(result, 1) << "Case 2: flip_performed expected for N atom on residue 11 chain A";
    }
    
    {
        int result = coot::pepflip_gemmi(st, "A", 9999, " ", "");
        EXPECT_EQ(result, 0) << "Case 3: flip_not_performed expected for non-existent residue 9999";
    }
    
    {
        int result = coot::pepflip_gemmi(st, "B", 10, " ", "");
        EXPECT_EQ(result, 0) << "Case 4: flip_not_performed expected for chain B (no residues there)";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}