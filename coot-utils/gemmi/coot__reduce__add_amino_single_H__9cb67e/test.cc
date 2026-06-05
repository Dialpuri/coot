#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static inline std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

static bool find_atom_by_name(const gemmi::Residue& res, const std::string& name, gemmi::Atom const** out) {
    std::string trimmed_name = trim(name);
    for (const auto& atom : res.atoms) {
        if (trim(atom.name) == trimmed_name) {
            *out = &atom;
            return true;
        }
    }
    return false;
}

TEST(OracleTestGemmi, add_amino_single_H) {
    // Load PDB into a gemmi::Structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    gemmi::Model* model = &st.models[0];
    ASSERT_NE(model, nullptr);

    ASSERT_FALSE(model->chains.empty());
    gemmi::Chain* chain = &model->chains[0];
    ASSERT_NE(chain, nullptr);
    ASSERT_FALSE(chain->residues.empty());

    // Case 1: Add H on backbone N using N, CA, C atoms
    {
        gemmi::Residue& residue = chain->residues[0]; // MET residue 1

        std::string H_at_name = "H1";
        std::string at_name_1 = "N";
        std::string at_name_2 = "CA";
        std::string at_name_3 = "C";
        double bl = 1.01;

        int atoms_before = static_cast<int>(residue.atoms.size());
        EXPECT_EQ(atoms_before, 19);

        coot::reduce::add_amino_single_H_gemmi(H_at_name, at_name_1, at_name_2, at_name_3, bl, residue);

        int atoms_after = static_cast<int>(residue.atoms.size());
        EXPECT_EQ(atoms_after, 20);

        // Check if H atom was added
        gemmi::Atom const* h_atom = nullptr;
        EXPECT_TRUE(find_atom_by_name(residue, "H1", &h_atom));
        if (h_atom) {
            EXPECT_NEAR(h_atom->pos.x, 13.6986, 1e-4);
            EXPECT_NEAR(h_atom->pos.y, 4.10144, 1e-4);
            EXPECT_NEAR(h_atom->pos.z, 73.7929, 1e-4);
        }
    }

    // Case 2: Edge case — non-existent atom name (XXX doesn't exist)
    {
        gemmi::Residue& residue2 = chain->residues[1]; // GLU residue 2
        ASSERT_NE(&residue2, nullptr);

        std::string H_at_name2 = "H99";
        std::string at_name_1_2 = "N";
        std::string at_name_2_2 = "CA";
        std::string at_name_3_2 = "XXX"; // non-existent atom
        double bl2 = 1.01;

        int atoms_before2 = static_cast<int>(residue2.atoms.size());
        EXPECT_EQ(atoms_before2, 15);

        coot::reduce::add_amino_single_H_gemmi(H_at_name2, at_name_1_2, at_name_2_2, at_name_3_2, bl2, residue2);

        int atoms_after2 = static_cast<int>(residue2.atoms.size());
        EXPECT_EQ(atoms_after2, 15);

        gemmi::Atom const* h_atom2 = nullptr;
        EXPECT_FALSE(find_atom_by_name(residue2, "H99", &h_atom2));
    }

    // Case 3: Additional coverage — add H on a side chain atom (e.g., OG on SER)
    // Using existing N, CA, CB atoms to add H on SER residue
    {
        gemmi::Residue& residue3 = chain->residues[34]; // SER residue 46
        ASSERT_NE(&residue3, nullptr);

        std::string H_at_name3 = "HG";
        std::string at_name_1_3 = "N";
        std::string at_name_2_3 = "CA";
        std::string at_name_3_3 = "CB";
        double bl3 = 1.01;

        int atoms_before3 = static_cast<int>(residue3.atoms.size());
        EXPECT_FALSE(atoms_before3 <= 0);

        coot::reduce::add_amino_single_H_gemmi(H_at_name3, at_name_1_3, at_name_2_3, at_name_3_3, bl3, residue3);

        int atoms_after3 = static_cast<int>(residue3.atoms.size());
        // With valid atoms, the hydrogen should be added
        EXPECT_GT(atoms_after3, atoms_before3);

        gemmi::Atom const* h_atom3 = nullptr;
        EXPECT_TRUE(find_atom_by_name(residue3, "HG", &h_atom3));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}