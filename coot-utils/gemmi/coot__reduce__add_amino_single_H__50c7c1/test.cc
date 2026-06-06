#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

static inline std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

TEST(OracleTestGemmi, add_amino_single_H) {
    // --- Case 1: Valid amino single H on backbone N of residue A/2 (GLU) ---
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        gemmi::Model& model = st.models[0];
        ASSERT_FALSE(model.chains.empty());

        gemmi::Chain* chain = nullptr;
        for (auto& ch : model.chains) {
            if (ch.name == "A") {
                chain = &ch;
                break;
            }
        }
        ASSERT_NE(chain, nullptr);

        // Get second residue (index 1, 0-indexed), A/2 GLU
        ASSERT_GE(static_cast<int>(chain->residues.size()), 2);
        gemmi::Residue& res_p = chain->residues[1];

        std::string H_at_name = "H";
        std::string first_neighb = "N";
        std::vector<std::string> second_neighb_vec = {"CA", "C"};
        double bl = 1.01;

        int atoms_before = static_cast<int>(res_p.atoms.size());
        EXPECT_EQ(atoms_before, 9);

        coot::reduce::add_amino_single_H_gemmi(H_at_name, first_neighb, second_neighb_vec, bl, res_p);

        int atoms_after = static_cast<int>(res_p.atoms.size());
        EXPECT_EQ(atoms_after, 10);
        EXPECT_TRUE(atoms_after > atoms_before);

        // Check if the H atom was actually created
        bool h_found = false;
        for (const auto& a : res_p.atoms) {
            if (trim(a.name) == trim(H_at_name)) {
                h_found = true;
                break;
            }
        }
        EXPECT_TRUE(h_found);
    }

    // --- Case 2: Invalid - second_neighb_vec size != 2, should do nothing ---
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        gemmi::Model& model = st.models[0];
        ASSERT_FALSE(model.chains.empty());

        gemmi::Chain* chain = nullptr;
        for (auto& ch : model.chains) {
            if (ch.name == "A") {
                chain = &ch;
                break;
            }
        }
        ASSERT_NE(chain, nullptr);

        // First residue A/1 MET
        ASSERT_GE(static_cast<int>(chain->residues.size()), 1);
        gemmi::Residue& res_p2 = chain->residues[0];

        std::string H_at_name = "H";
        std::string first_neighb = "N";
        std::vector<std::string> bad_second_neighb_vec = {"CA"}; // size 1, not 2
        double bl = 1.01;

        int atoms_before2 = static_cast<int>(res_p2.atoms.size());
        EXPECT_EQ(atoms_before2, 8);

        coot::reduce::add_amino_single_H_gemmi(H_at_name, first_neighb, bad_second_neighb_vec, bl, res_p2);

        int atoms_after2 = static_cast<int>(res_p2.atoms.size());
        EXPECT_EQ(atoms_after2, 8);
        EXPECT_FALSE(atoms_after2 > atoms_before2);
    }

    // --- Case 3: Empty second_neighb_vec (size 0), should also do nothing ---
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        gemmi::Model& model = st.models[0];
        ASSERT_FALSE(model.chains.empty());

        gemmi::Chain* chain = nullptr;
        for (auto& ch : model.chains) {
            if (ch.name == "A") {
                chain = &ch;
                break;
            }
        }
        ASSERT_NE(chain, nullptr);

        // A/1 MET
        ASSERT_GE(static_cast<int>(chain->residues.size()), 1);
        gemmi::Residue& res_p3 = chain->residues[0];

        std::string H_at_name = "H";
        std::string first_neighb = "N";
        std::vector<std::string> empty_second_neighb_vec; // size 0
        double bl = 1.01;

        int atoms_before3 = static_cast<int>(res_p3.atoms.size());
        EXPECT_EQ(atoms_before3, 8);

        coot::reduce::add_amino_single_H_gemmi(H_at_name, first_neighb, empty_second_neighb_vec, bl, res_p3);

        int atoms_after3 = static_cast<int>(res_p3.atoms.size());
        EXPECT_EQ(atoms_after3, 8);
        EXPECT_FALSE(atoms_after3 > atoms_before3);
    }
}