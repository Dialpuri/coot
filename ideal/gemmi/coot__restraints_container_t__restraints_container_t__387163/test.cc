#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, restraints_container_constructor) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];
    
    std::vector<std::pair<bool, gemmi::Residue*>> residues;
    int n_to_take = std::min((int)chain.residues.size(), 3);
    
    for (int i = 0; i < n_to_take; i++) {
        gemmi::Residue* res = &chain.residues[i];
        if (res) {
            residues.push_back(std::make_pair(false, res));
        }
    }
    
    clipper::Xmap<float> dummy_xmap;
    
    coot::restraints_container_t_gemmi restraints(residues, "test", &st, &dummy_xmap);
    
    EXPECT_EQ(restraints.get_residues_vec_size(), 3);
    EXPECT_TRUE(restraints.has_valid_residues());
}

TEST(OracleTest, restraints_container_empty) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    std::vector<std::pair<bool, gemmi::Residue*>> residues;
    
    clipper::Xmap<float> dummy_xmap;
    
    coot::restraints_container_t_gemmi restraints(residues, "test", &st, &dummy_xmap);
    
    EXPECT_EQ(restraints.get_residues_vec_size(), 0);
    EXPECT_FALSE(restraints.has_valid_residues());
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}