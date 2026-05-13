#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, NextResidueNumberInChain) {
    std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
    gemmi::Structure st = gemmi::read_pdb_file(pdb_path);
    
    ASSERT_FALSE(st.models.empty());
    
    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());
    
    gemmi::Chain& chain = model.chains[0];
    std::cout << "CHAIN ID: " << chain.name << std::endl;
    std::cout << "NUMBER OF RESIDUES: " << chain.residues.size() << std::endl;
    
    size_t nres = chain.residues.size();
    ASSERT_GT(nres, 0);
    
    const gemmi::Residue* first_res = nullptr;
    const gemmi::Residue* last_res = nullptr;
    int min_seq = 99999;
    int max_seq = -99999;
    
    for (const gemmi::Residue& res : chain.residues) {
        int seq = res.seqid.num.value;
        if (seq < min_seq) { min_seq = seq; first_res = &res; }
        if (seq > max_seq) { max_seq = seq; last_res = &res; }
    }
    
    std::cout << "FIRST RESidue SEQNUM: " << (first_res ? first_res->seqid.num.value : -1) << std::endl;
    std::cout << "LAST RESidue SEQNUM: " << (last_res ? last_res->seqid.num.value : -1) << std::endl;
    
    bool new_res_no_by_hundreds = false;
    std::pair<short int, int> result1 = coot::molecule_t::next_residue_number_in_chain_gemmi(chain, new_res_no_by_hundreds);
    std::cout << "INPUT new_res_no_by_hundreds: " << new_res_no_by_hundreds << std::endl;
    std::cout << "OUTPUT result.first (state): " << result1.first << std::endl;
    std::cout << "OUTPUT result.second (next_resno): " << result1.second << std::endl;
    
    new_res_no_by_hundreds = true;
    std::pair<short int, int> result2 = coot::molecule_t::next_residue_number_in_chain_gemmi(chain, new_res_no_by_hundreds);
    std::cout << "INPUT new_res_no_by_hundreds: " << new_res_no_by_hundreds << std::endl;
    std::cout << "OUTPUT result.first (state): " << result2.first << std::endl;
    std::cout << "OUTPUT result.second (next_resno): " << result2.second << std::endl;
    
    gemmi::Chain empty_chain;
    empty_chain.name = " ";
    std::pair<short int, int> result3 = coot::molecule_t::next_residue_number_in_chain_gemmi(empty_chain, false);
    std::cout << "INPUT chain (empty): true" << std::endl;
    std::cout << "OUTPUT result.first (state): " << result3.first << std::endl;
    std::cout << "OUTPUT result.second (next_resno): " << result3.second << std::endl;
    
    EXPECT_EQ(result1.first, 1);
    EXPECT_EQ(result1.second, 299);
    
    EXPECT_EQ(result2.first, 1);
    EXPECT_EQ(result2.second, 301);
    
    EXPECT_EQ(result3.first, 0);
    EXPECT_EQ(result3.second, 1);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}