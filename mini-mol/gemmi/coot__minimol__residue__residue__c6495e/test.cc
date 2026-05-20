#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ResidueConstructorWithFilters) {
    // Load a PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_EQ(st.models.size(), 1) << "Failed to read PDB file";

    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];
    gemmi::Residue& residue = chain.residues[0];

    // Test 1: Constructor with empty filter vector (should include all atoms)
    {
        std::vector<std::string> keep_only_these_atoms;
        
        coot::minimol::residue r(&residue, keep_only_these_atoms);
        EXPECT_EQ(r.seqnum, 1);
        EXPECT_EQ(r.ins_code, "");
        EXPECT_EQ(r.name, "MET");
        EXPECT_EQ(r.n_atoms(), 0);
    }

    // Test 2: Constructor with filter - keep only CA atoms
    {
        std::vector<std::string> keep_only_these_atoms;
        keep_only_these_atoms.push_back(" CA ");
        
        coot::minimol::residue r(&residue, keep_only_these_atoms);
        EXPECT_EQ(r.seqnum, 1);
        EXPECT_EQ(r.ins_code, "");
        EXPECT_EQ(r.name, "MET");
        EXPECT_EQ(r.n_atoms(), 1);
        
        if (r.n_atoms() > 0) {
            EXPECT_EQ(r[0].name, " CA ");
        }
    }

    // Test 3: Constructor with filter that doesn't match any atoms
    {
        std::vector<std::string> keep_only_these_atoms;
        keep_only_these_atoms.push_back(" ZZ ");
        
        coot::minimol::residue r(&residue, keep_only_these_atoms);
        EXPECT_EQ(r.seqnum, 1);
        EXPECT_EQ(r.n_atoms(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}