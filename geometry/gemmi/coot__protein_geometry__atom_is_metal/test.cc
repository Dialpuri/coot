#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, AtomIsMetal) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get the first model
    ASSERT_FALSE(st.models.empty());
    const gemmi::Model& model = st.models[0];
    
    // Get first chain
    ASSERT_FALSE(model.chains.empty());
    const gemmi::Chain& chain = model.chains[0];
    
    // Get first residue
    ASSERT_FALSE(chain.residues.empty());
    const gemmi::Residue& res = chain.residues[0];
    
    // Get atoms from the residue
    const std::vector<gemmi::Atom>& atoms = res.atoms;
    ASSERT_FALSE(atoms.empty()) << "No atoms in first residue";
    
    // Test case 1: Test with first atom (N) - should NOT be metal
    {
        const gemmi::Atom& atom = atoms[0];
        std::string atom_name = atom.name;
        bool is_metal = coot::protein_geometry::atom_is_metal_gemmi(atom);
        std::cout << "INPUT  atom_name: " << atom_name << std::endl;
        std::cout << "OUTPUT is_metal: " << (is_metal ? "true" : "false") << std::endl;
        EXPECT_FALSE(is_metal) << "Atom " << atom_name << " should not be a metal";
    }
    
    // Test case 2: Check second atom in the same residue (CA)
    if (atoms.size() > 1) {
        const gemmi::Atom& atom = atoms[1];
        std::string atom_name = atom.name;
        bool is_metal = coot::protein_geometry::atom_is_metal_gemmi(atom);
        std::cout << "INPUT  atom_name: " << atom_name << std::endl;
        std::cout << "OUTPUT is_metal: " << (is_metal ? "true" : "false") << std::endl;
        EXPECT_FALSE(is_metal) << "Atom " << atom_name << " should not be a metal";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}