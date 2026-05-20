#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, AddHydrogenAtom) {
    // Load the test PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_FALSE(st.models.empty());
    
    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());
    
    gemmi::Chain& chain = model.chains[0];
    ASSERT_FALSE(chain.residues.empty());
    
    gemmi::Residue& res = chain.residues[0];
    
    // Check existing atoms
    int n_atoms_before = res.atoms.size();
    EXPECT_EQ(n_atoms_before, 19);
    EXPECT_EQ(res.name, "MET");
    EXPECT_EQ(res.seqid.num.value, 1);

    // Case 1: Add a new hydrogen atom
    {
        clipper::Coord_orth pos(1.0, 2.0, 3.0);
        float bf = 20.0f;
        std::string altconf = "";
        
        gemmi::Atom* new_H = coot::reduce::add_hydrogen_atom_gemmi(" HA ", pos, bf, altconf, res);
        
        EXPECT_TRUE(new_H != nullptr);
        if (new_H) {
            EXPECT_EQ(new_H->name, std::string(" HA "));
            EXPECT_NEAR(new_H->pos.x, 1.0, 1e-4);
            EXPECT_NEAR(new_H->pos.y, 2.0, 1e-4);
            EXPECT_NEAR(new_H->pos.z, 3.0, 1e-4);
        }
    }

    // Case 2: Add another hydrogen with same name (should update existing)
    {
        clipper::Coord_orth pos(1.5, 2.5, 3.5);
        float bf = 25.0f;
        std::string altconf = "";
        
        gemmi::Atom* existing_H = coot::reduce::add_hydrogen_atom_gemmi(" HA ", pos, bf, altconf, res);
        
        EXPECT_TRUE(existing_H != nullptr);
        if (existing_H) {
            EXPECT_EQ(existing_H->name, std::string(" HA "));
            EXPECT_NEAR(existing_H->pos.x, 1.5, 1e-4);
            EXPECT_NEAR(existing_H->pos.y, 2.5, 1e-4);
            EXPECT_NEAR(existing_H->pos.z, 3.5, 1e-4);
        }
    }

    // Check final atom count - should still be 19 since we updated existing
    int n_atoms_after = res.atoms.size();
    EXPECT_EQ(n_atoms_after, 19);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}