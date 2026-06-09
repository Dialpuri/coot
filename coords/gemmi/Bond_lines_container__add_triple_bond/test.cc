#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "function.hh"

TEST(OracleTest, add_triple_bond) {
    // Read the PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get model 0 (gemmi is 0-indexed; MMDB model 1 = gemmi model 0)
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Get chain A
    const gemmi::Chain* chain = model.find_chain("A");
    ASSERT_NE(chain, nullptr);

    // Get first residue
    ASSERT_FALSE(chain->residues.empty());
    const gemmi::Residue& res = chain->residues[0];

    EXPECT_EQ(chain->name, "A");
    EXPECT_EQ(res.seqid.num.value, 1);
    EXPECT_EQ(res.name, "MET");

    int n_atoms = static_cast<int>(res.atoms.size());
    EXPECT_EQ(n_atoms, 8);

    // Build atom vector for the residue
    std::vector<gemmi::Atom> atoms;
    for (auto& atom : res.atoms) {
        atoms.push_back(atom);
    }

    // Set up bond storage (mirrors Bond_lines_container::bonds)
    std::map<int, Bond_lines> bonds_storage;

    // Setup parameters
    int imol = 0;
    int imodel = 1;
    int iat_1 = 0;
    int iat_2 = (n_atoms > 1) ? 1 : 0;
    int atom_colour_type = 0;
    int udd_user_defined_atom_colour_index_handle = -1;

    EXPECT_EQ(iat_1, 0);
    EXPECT_EQ(iat_2, 1);
    EXPECT_EQ(imol, 0);
    EXPECT_EQ(imodel, 1);

    // Count total bonds before
    int total_bonds_before = 0;
    for (auto& kv : bonds_storage) {
        total_bonds_before += kv.second.size();
    }

    // Create protein geometry for neighbour normal
    coot::protein_geometry geom;
    geom.init_standard();

    // Call the gemmi port
    add_triple_bond_gemmi(
        bonds_storage,
        imol,
        imodel,
        iat_1,
        iat_2,
        atoms,
        res.name,
        atom_colour_type,
        nullptr,
        udd_user_defined_atom_colour_index_handle,
        &geom);

    // Count total bonds after
    int total_bonds_after = 0;
    for (auto& kv : bonds_storage) {
        total_bonds_after += kv.second.size();
    }

    // A triple bond adds 6 bond segments
    EXPECT_EQ(total_bonds_after - total_bonds_before, 6);

    // Case 2: Check element comparison with gemmi
    // In MMDB the original code compared char* pointers (broken).
    // In gemmi, element.name() returns std::string, so comparison is by content.
    int c1 = -1, c2 = -1;
    for (int i = 0; i < n_atoms && c2 < 0; i++) {
        for (int j = i + 1; j < n_atoms && c2 < 0; j++) {
            if (atoms[i].element.name() == atoms[j].element.name()) {
                c1 = i;
                c2 = j;
            }
        }
    }

    // With gemmi string comparison, we should find atoms with the same element
    // MET has multiple C and H atoms, so this should succeed
    EXPECT_TRUE(c1 >= 0 && c2 >= 0);

    // Verify the bonds were added to the correct colour buckets
    EXPECT_FALSE(bonds_storage.empty());
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
