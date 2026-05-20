#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, MoleculeDiameter) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    std::vector<gemmi::Atom> atoms;
    for (auto& model : st.models)
        for (auto& chain : model.chains)
            for (auto& res : chain.residues)
                for (auto& atom : res.atoms)
                    atoms.push_back(atom);
    
    float diameter = coot::get_molecule_diameter_gemmi(atoms);
    EXPECT_NEAR(diameter, 51.6117, 1e-4);
    
    EXPECT_EQ(atoms.size(), 2107);
}

// Additional test case to cover alternative branch (s.size() <= 10)
// With 1 atom, no valid pairs exist (idx_1 != idx_2 is never true)
TEST(MoleculeDiameterSmallSet, SmallDiameter) {
    std::vector<gemmi::Atom> small_atoms;
    gemmi::Atom at;
    at.pos = gemmi::Position(0, 0, 0);
    small_atoms.push_back(at);
    
    float diameter = coot::get_molecule_diameter_gemmi(small_atoms);
    // With only 1 atom, no valid pairs, s.size() = 0 <= 10, function returns -1
    EXPECT_EQ(diameter, -1.0f);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}