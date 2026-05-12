#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper to count TER atoms in a residue (by checking for "TER" in atom name)
static int count_TER_atoms(const gemmi::Residue& residue) {
    int count = 0;
    for (const auto& atom : residue.atoms) {
        if (atom.name.find("TER") != std::string::npos) {
            count++;
        }
    }
    return count;
}

TEST(OracleTest, RemoveTERInternal) {
    // case: setup - load the example PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Use first residue in chain A as test residue
    gemmi::Residue* res = nullptr;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name == "A" && !chain.residues.empty()) {
                res = &chain.residues[0];
                break;
            }
        }
        if (res) break;
    }
    ASSERT_NE(res, nullptr) << "Residue not found in chain A";
    
    // case: before - should have 0 TER atoms initially
    int ter_count_before = count_TER_atoms(*res);
    EXPECT_EQ(ter_count_before, 0);
    
    // Add a TER atom so remove_TER_internal has something to delete
    gemmi::Atom new_ter_atom;
    if (!res->atoms.empty()) {
        new_ter_atom = res->atoms[0];  // Copy position from first atom
        new_ter_atom.name = "TER ";
        // Add to residue
        res->atoms.push_back(new_ter_atom);
    }
    
    // case: after add - should have 1 TER atom
    int ter_count_after_add = count_TER_atoms(*res);
    EXPECT_EQ(ter_count_after_add, 1);
    
    // Now call remove_TER_internal
    coot::remove_TER_internal_gemmi(st, res);
    
    // case: after remove - should have 0 TER atoms
    int ter_count_after_remove = count_TER_atoms(*res);
    EXPECT_EQ(ter_count_after_remove, 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}