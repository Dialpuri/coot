#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, setHasHydrogenAtomsState) {
    // case: basic functionality
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
        ASSERT_FALSE(st.models.empty()) << "no models";
        
        // Get the first chain and residue
        const auto& model = st.models[0];
        ASSERT_FALSE(model.chains.empty()) << "no chains";
        const auto& chain = model.chains[0];
        ASSERT_FALSE(chain.residues.empty()) << "no residues";
        const auto& residue = chain.residues[0];
        
        // Test with a residue that should have hydrogens
        bool model_has_hydrogen_atoms = true; // initial state
        coot::set_has_hydrogen_atoms_state(residue.atoms, model_has_hydrogen_atoms);
        EXPECT_TRUE(model_has_hydrogen_atoms) << "should have hydrogens";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}