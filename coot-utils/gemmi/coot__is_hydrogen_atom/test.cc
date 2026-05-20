#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, IsHydrogenAtom) {
    const std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb";
    gemmi::Structure st = gemmi::read_pdb_file(pdb_path);
    
    // First, let's dump some atoms to see what names we have
    int atom_count = 0;
    for (auto& model : st.models)
        for (auto& chain : model.chains)
            for (auto& res : chain.residues) {
                for (const auto& atom : res.atoms) {
                    if (atom_count < 5) {
                        atom_count++;
                    }
                }
                if (atom_count >= 5) break;
            }
    
    // case 1: Check for hydrogen atoms by element
    {
        bool found_hydrogen = false;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    for (const auto& atom : res.atoms) {
                        if (atom.is_hydrogen()) {
                            bool result = coot::is_hydrogen_atom_gemmi(atom);
                            EXPECT_TRUE(result) << "is_hydrogen should be true for hydrogen atom " << atom.name;
                            found_hydrogen = true;
                            break;
                        }
                    }
        EXPECT_TRUE(found_hydrogen) << "Should find at least one hydrogen atom";
    }
    
    // case 2: Check for non-hydrogen atoms
    {
        bool found_non_hydrogen = false;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    for (const auto& atom : res.atoms) {
                        if (!atom.is_hydrogen()) {
                            bool result = coot::is_hydrogen_atom_gemmi(atom);
                            EXPECT_FALSE(result) << "is_hydrogen should be false for non-hydrogen atom " << atom.name;
                            found_non_hydrogen = true;
                            break;
                        }
                    }
        EXPECT_TRUE(found_non_hydrogen) << "Should find at least one non-hydrogen atom";
    }
    
    // case 3: count all hydrogen vs non-hydrogen atoms
    {
        int hydrogen_count = 0;
        int non_hydrogen_count = 0;
        
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    for (const auto& atom : res.atoms) {
                        bool result = coot::is_hydrogen_atom_gemmi(atom);
                        bool is_hydrogen_elem = atom.is_hydrogen();
                        
                        if (is_hydrogen_elem) {
                            EXPECT_TRUE(result) << "Atom " << atom.name 
                                                << " is hydrogen but is_hydrogen_atom_gemmi returned false";
                            hydrogen_count++;
                        } else {
                            EXPECT_FALSE(result) << "Atom " << atom.name 
                                                 << " is not hydrogen but is_hydrogen_atom_gemmi returned true";
                            non_hydrogen_count++;
                        }
                    }
        
        EXPECT_GT(hydrogen_count, 0) << "Should find at least one hydrogen atom";
        EXPECT_GT(non_hydrogen_count, 0) << "Should find at least one non-hydrogen atom";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}