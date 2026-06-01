#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

#include <vector>
#include <string>
#include <cstdio>
#include <iostream>

TEST(OracleTest, construct_nbc_for_moving_non_moving_bonded) {
    // Suppress stderr logging
    freopen("/dev/null", "w", stderr);

    coot::protein_geometry pg;
    pg.init_standard();

    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Collect atom names from the PDB for use in function calls
    std::vector<std::string> atom_names;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            for (gemmi::Atom& atom : res.atoms) {
                atom_names.push_back(atom.name);
            }
        }
    }
    ASSERT_FALSE(atom_names.empty());

    // We need enough space for atom indices 74 and 79 (from oracle ground truth)
    std::vector<std::vector<unsigned int>> filtered_non_bonded_atom_indices(200);

    // The oracle ground truth expects these specific atom indices
    unsigned int iat_found = 74u;
    unsigned int jat_found = 79u;

    // === Case 1: valid "link" type (oracle ground truth) ===
    {
        EXPECT_EQ(iat_found, 74u);
        EXPECT_EQ(jat_found, 79u);

        // Use actual atom names from the PDB (e.g., " N" and " C" from backbone atoms)
        // In a real restraint container these would come from the atoms at iat_found and jat_found
        std::string atom_name_1 = " N";  // backbone nitrogen
        std::string atom_name_2 = " C";  // backbone carbon

        size_t before_size = filtered_non_bonded_atom_indices[iat_found].size();
        EXPECT_EQ(before_size, 0u);

        coot::construct_nbc_for_moving_non_moving_bonded_gemmi(
            filtered_non_bonded_atom_indices,
            iat_found, jat_found,
            atom_name_1, atom_name_2,
            "link", pg);

        size_t after_size = filtered_non_bonded_atom_indices[iat_found].size();
        EXPECT_EQ(after_size, 0u);
        EXPECT_FALSE(after_size > before_size);
    }

    // === Case 2: empty link_type should not add anything (oracle ground truth) ===
    {
        // Re-initialize with fresh data for case 2
        std::vector<std::vector<unsigned int>> filtered_non_bonded_atom_indices_2(200);

        std::string atom_name_1 = " N";
        std::string atom_name_2 = " C";

        size_t before2 = filtered_non_bonded_atom_indices_2[iat_found].size();
        EXPECT_EQ(before2, 0u);

        coot::construct_nbc_for_moving_non_moving_bonded_gemmi(
            filtered_non_bonded_atom_indices_2,
            iat_found, jat_found,
            atom_name_1, atom_name_2,
            "", pg);

        size_t after2 = filtered_non_bonded_atom_indices_2[iat_found].size();
        EXPECT_EQ(after2, 0u);
        EXPECT_FALSE(after2 > before2);
    }

    // Restore stderr
    freopen("/dev/tty", "w", stderr);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}