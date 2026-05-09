#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(RemoveTERInternal, ChainA) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    auto& model = st.models[0];
    gemmi::Chain* chain_a = nullptr;
    for (auto& c : model.chains) {
        if (c.name == "A") {
            chain_a = &c;
            break;
        }
    }
    ASSERT_NE(chain_a, nullptr);
    
    int n_res_before = chain_a->residues.size();
    EXPECT_EQ(n_res_before, 267);
    
    auto& last_res = chain_a->residues.back();
    
    coot::molecule_t::remove_TER_internal_gemmi(last_res);
    
    int n_atoms = last_res.atoms.size();
    EXPECT_EQ(n_atoms, 9);
}

TEST(RemoveTERInternal, ChainB) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    auto& model = st.models[0];
    gemmi::Chain* chain_b = nullptr;
    for (auto& c : model.chains) {
        if (c.name == "B") {
            chain_b = &c;
            break;
        }
    }
    ASSERT_NE(chain_b, nullptr);
    
    int n_res_before = chain_b->residues.size();
    EXPECT_EQ(n_res_before, 1);
    
    auto& last_res = chain_b->residues.back();
    
    coot::molecule_t::remove_TER_internal_gemmi(last_res);
    
    int n_atoms = last_res.atoms.size();
    EXPECT_EQ(n_atoms, 1);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}