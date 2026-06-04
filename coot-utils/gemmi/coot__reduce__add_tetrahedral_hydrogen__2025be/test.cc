#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper to find a chain by name
static gemmi::Chain* find_chain(gemmi::Model& model, const std::string& chain_name) {
    for (auto& chain : model.chains) {
        if (chain.name == chain_name) return &chain;
    }
    return nullptr;
}

// Helper to find a residue by its index in the chain's residue vector
static gemmi::Residue* find_residue_by_index(gemmi::Chain& chain, int index) {
    if (index < 0 || index >= static_cast<int>(chain.residues.size())) return nullptr;
    return &chain.residues[index];
}

// Helper to find an atom by name (trimmed) and altconf in a residue
static const gemmi::Atom* find_atom(const gemmi::Residue& residue, const std::string& atom_name, const std::string& altconf) {
    for (const auto& atom : residue.atoms) {
        // altconf "" means default: match ' ' or '\0'
        if (!altconf.empty() && atom.altloc != altconf[0]) continue;
        if (altconf.empty() && atom.altloc != ' ' && atom.altloc != '\0') continue;
        // Trim the atom name for comparison
        std::string trimmed_name = atom.name;
        size_t first = trimmed_name.find_first_not_of(" \t\n\r");
        size_t last = trimmed_name.find_last_not_of(" \t\n\r");
        if (first == std::string::npos) continue;
        trimmed_name = trimmed_name.substr(first, last - first + 1);
        std::string target = atom_name;
        size_t first_t = target.find_first_not_of(" \t\n\r");
        size_t last_t = target.find_last_not_of(" \t\n\r");
        if (first_t != std::string::npos)
            target = target.substr(first_t, last_t - first_t + 1);
        if (trimmed_name == target) return &atom;
    }
    return nullptr;
}

TEST(OracleTest, add_tetrahedral_hydrogen) {
    // Load the PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Get model 1 (0-indexed in gemmi)
    gemmi::Model& model = st.models[0];

    // Get chain A
    gemmi::Chain* chain = find_chain(model, "A");
    ASSERT_NE(chain, nullptr);

    // Case 1: VAL residue at index 6 (seq 7) - has CA, CB, CG1, CG2
    // Use unique H name "HBX" so atoms_before != atoms_after
    {
        gemmi::Residue* res = find_residue_by_index(*chain, 6);
        ASSERT_NE(res, nullptr);

        EXPECT_STREQ(res->name.c_str(), "VAL");

        int num_atoms_before = static_cast<int>(res->atoms.size());
        EXPECT_EQ(num_atoms_before, 16);

        std::string H_at_name = "HBX";
        std::string at_central_name = "CB";
        std::string neighb_at_name_1 = "CA";
        std::string neighb_at_name_2 = "CG1";
        std::string neighb_at_name_3 = "CG2";
        double bond_length = 1.0;

        coot::reduce::add_tetrahedral_hydrogen_gemmi(H_at_name, at_central_name,
                                                     neighb_at_name_1, neighb_at_name_2, neighb_at_name_3,
                                                     bond_length, *res);

        int num_atoms_after = static_cast<int>(res->atoms.size());
        EXPECT_EQ(num_atoms_after, 17);

        const gemmi::Atom* hb_atom = find_atom(*res, "HBX", "");
        EXPECT_NE(hb_atom, nullptr);
        if (hb_atom) {
            EXPECT_NEAR(hb_atom->pos.x, 26.791, 1e-4);
            EXPECT_NEAR(hb_atom->pos.y, 8.70356, 1e-4);
            EXPECT_NEAR(hb_atom->pos.z, 75.2085, 1e-4);
        }
    }

    // Case 2: Invalid neighbor (NONEXISTENT) - should NOT add H
    {
        gemmi::Residue* res = find_residue_by_index(*chain, 6);
        ASSERT_NE(res, nullptr);

        EXPECT_STREQ(res->name.c_str(), "VAL");

        int num_atoms_before = static_cast<int>(res->atoms.size());
        EXPECT_EQ(num_atoms_before, 17);

        std::string H_at_name = "HBX2";
        std::string at_central_name = "CB";
        std::string neighb_at_name_1 = "CA";
        std::string neighb_at_name_2 = "CG1";
        std::string neighb_at_name_3 = "NONEXISTENT";
        double bond_length = 1.0;

        coot::reduce::add_tetrahedral_hydrogen_gemmi(H_at_name, at_central_name,
                                                     neighb_at_name_1, neighb_at_name_2, neighb_at_name_3,
                                                     bond_length, *res);

        int num_atoms_after = static_cast<int>(res->atoms.size());
        EXPECT_EQ(num_atoms_after, 17);

        const gemmi::Atom* hb_atom = find_atom(*res, "HBX2", "");
        EXPECT_EQ(hb_atom, nullptr);
    }

    // Case 3: THR residue at index 13 (seq 14) - has CB, CA, OG1, CG2
    // Use unique H name "HBX3"
    {
        gemmi::Residue* res = find_residue_by_index(*chain, 13);
        ASSERT_NE(res, nullptr);

        EXPECT_STREQ(res->name.c_str(), "THR");

        int num_atoms_before = static_cast<int>(res->atoms.size());
        EXPECT_EQ(num_atoms_before, 14);

        std::string H_at_name = "HBX3";
        std::string at_central_name = "CB";
        std::string neighb_at_name_1 = "CA";
        std::string neighb_at_name_2 = "OG1";
        std::string neighb_at_name_3 = "CG2";
        double bond_length = 1.0;

        coot::reduce::add_tetrahedral_hydrogen_gemmi(H_at_name, at_central_name,
                                                     neighb_at_name_1, neighb_at_name_2, neighb_at_name_3,
                                                     bond_length, *res);

        int num_atoms_after = static_cast<int>(res->atoms.size());
        EXPECT_EQ(num_atoms_after, 15);

        const gemmi::Atom* hb_atom = find_atom(*res, "HBX3", "");
        EXPECT_NE(hb_atom, nullptr);
        if (hb_atom) {
            EXPECT_NEAR(hb_atom->pos.x, 17.4084, 1e-4);
            EXPECT_NEAR(hb_atom->pos.y, 13.9789, 1e-4);
            EXPECT_NEAR(hb_atom->pos.z, 55.2145, 1e-4);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}