#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a residue pointer in a model by chain name and residue sequence number
// Returns {Residue*, Chain*} pair
static std::pair<gemmi::Residue*, gemmi::Chain*> find_residue(gemmi::Model& model, const std::string& chain_name, int seqnum) {
    for (auto& chain : model.chains) {
        if (chain.name == chain_name) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == seqnum) {
                    return {&res, &chain};
                }
            }
        }
    }
    return {nullptr, nullptr};
}

// Helper: find an atom pointer in a residue by name (trimmed comparison)
static gemmi::Atom* find_atom(gemmi::Residue& res, const std::string& name) {
    std::string trimmed = name;
    // Trim whitespace from both ends
    size_t first = trimmed.find_first_not_of(" \t");
    if (first == std::string::npos) return nullptr;
    size_t last = trimmed.find_last_not_of(" \t");
    trimmed = trimmed.substr(first, last - first + 1);

    for (auto& atom : res.atoms) {
        std::string atom_trimmed = atom.name;
        size_t af = atom_trimmed.find_first_not_of(" \t");
        if (af == std::string::npos) continue;
        size_t al = atom_trimmed.find_last_not_of(" \t");
        atom_trimmed = atom_trimmed.substr(af, al - af + 1);
        if (atom_trimmed == trimmed) return &atom;
    }
    return nullptr;
}

TEST(OracleTest, add_amino_hydrogens) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model& model = st.models[0];

    // ---- Case 1: ASN residue 3 — 6-param version with valid neighbors (hydrogens already exist) ----
    {
        auto [residue_p, chain_p] = find_residue(model, "A", 3);
        ASSERT_NE(residue_p, nullptr);
        ASSERT_NE(chain_p, nullptr);
        EXPECT_EQ(chain_p->name[0], 'A');
        EXPECT_EQ(residue_p->seqid.num.value, 3);
        EXPECT_EQ(residue_p->name, "ASN");

        int atoms_before = static_cast<int>(residue_p->atoms.size());
        EXPECT_EQ(atoms_before, 14);

        gemmi::Atom* h1 = find_atom(*residue_p, "HD21");
        gemmi::Atom* h2 = find_atom(*residue_p, "HD22");
        EXPECT_NE(h1, nullptr);
        EXPECT_NE(h2, nullptr);

        std::string H_at_name_1 = "HD21";
        std::string H_at_name_2 = "HD22";
        std::string first_neighb = "ND2";
        std::vector<std::string> second_neighb_vec;
        second_neighb_vec.push_back("CG");
        std::map<std::string, std::vector<std::string>> third_neighb_map;
        std::vector<std::string> thirds;
        thirds.push_back("OD1");
        third_neighb_map["CG"] = thirds;
        double bl_amino = 1.01;

        EXPECT_NO_THROW(
            coot::reduce::add_amino_hydrogens_gemmi(H_at_name_1, H_at_name_2, first_neighb,
                                        second_neighb_vec, third_neighb_map, bl_amino, *residue_p, st, chain_p)
        );

        int atoms_after = static_cast<int>(residue_p->atoms.size());
        EXPECT_EQ(atoms_after, 14);
    }

    // ---- Case 2: GLN residue 5 ----
    {
        auto [residue_p, chain_p] = find_residue(model, "A", 5);
        ASSERT_NE(residue_p, nullptr);
        ASSERT_NE(chain_p, nullptr);

        EXPECT_EQ(residue_p->seqid.num.value, 5);
        EXPECT_EQ(residue_p->name, "GLN");

        int atoms_before = static_cast<int>(residue_p->atoms.size());
        EXPECT_EQ(atoms_before, 17);

        gemmi::Atom* he21 = find_atom(*residue_p, "HE21");
        EXPECT_NE(he21, nullptr);

        std::string H_at_name_1 = "HE21";
        std::string H_at_name_2 = "HE22";
        std::string first_neighb = "NE2";
        std::vector<std::string> second_neighb_vec;
        second_neighb_vec.push_back("CD");
        std::map<std::string, std::vector<std::string>> third_neighb_map;
        std::vector<std::string> thirds;
        thirds.push_back("OE1");
        third_neighb_map["CD"] = thirds;
        double bl_amino = 1.01;

        EXPECT_NO_THROW(
            coot::reduce::add_amino_hydrogens_gemmi(H_at_name_1, H_at_name_2, first_neighb,
                                        second_neighb_vec, third_neighb_map, bl_amino, *residue_p, st, chain_p)
        );

        int atoms_after = static_cast<int>(residue_p->atoms.size());
        EXPECT_EQ(atoms_after, 17);
    }

    // ---- Case 3: edge case — empty second_neighb_vec (no-op) ----
    {
        auto [residue_p, chain_p] = find_residue(model, "A", 11);
        ASSERT_NE(residue_p, nullptr);
        ASSERT_NE(chain_p, nullptr);

        EXPECT_EQ(residue_p->seqid.num.value, 11);
        EXPECT_EQ(residue_p->name, "GLY");

        int atoms_before = static_cast<int>(residue_p->atoms.size());
        EXPECT_EQ(atoms_before, 7);

        std::vector<std::string> second_neighb_vec;
        std::map<std::string, std::vector<std::string>> third_neighb_map;

        EXPECT_NO_THROW(
            coot::reduce::add_amino_hydrogens_gemmi("HXXX", "HYYY", "N",
                                        second_neighb_vec, third_neighb_map, 1.01, *residue_p, st, chain_p)
        );

        int atoms_after = static_cast<int>(residue_p->atoms.size());
        EXPECT_EQ(atoms_after, 7);
    }

    // ---- Case 4: 6-param overload with empty second_neighb_vec on a different residue ----
    // Tests that the function handles the case where second_neighb_vec is empty but
    // first_neighb is valid - should be a no-op since second neighbors are required.
    {
        auto [residue_p, chain_p] = find_residue(model, "A", 7);
        // residue 7 might or might not exist - check
        if (residue_p) {
            int atoms_before = static_cast<int>(residue_p->atoms.size());

            std::vector<std::string> second_neighb_vec;
            std::map<std::string, std::vector<std::string>> third_neighb_map;

            EXPECT_NO_THROW(
                coot::reduce::add_amino_hydrogens_gemmi("HN", "HXXX", "N",
                                            second_neighb_vec, third_neighb_map, 1.01, *residue_p, st, chain_p)
            );

            int atoms_after = static_cast<int>(residue_p->atoms.size());
            EXPECT_EQ(atoms_after, atoms_before);
        }
    }

    // ---- Case 5: 6-param overload with non-empty second_neighb but nonexistent atom names ----
    // Tests with a residue where the target hydrogen atoms don't exist but neighbors do.
    // The function should still succeed without crashing even if no hydrogens are added.
    {
        auto [residue_p, chain_p] = find_residue(model, "A", 3);
        ASSERT_NE(residue_p, nullptr);
        ASSERT_NE(chain_p, nullptr);

        int atoms_before = static_cast<int>(residue_p->atoms.size());

        std::string H_at_name_1 = "HZ1";
        std::string H_at_name_2 = "HZ2";
        std::string first_neighb = "NZ";  // ASN doesn't have NZ
        std::vector<std::string> second_neighb_vec;
        second_neighb_vec.push_back("CZ"); // ASN doesn't have CZ either
        std::map<std::string, std::vector<std::string>> third_neighb_map;
        std::vector<std::string> thirds;
        thirds.push_back("OZ1");
        third_neighb_map["CZ"] = thirds;
        double bl_amino = 1.01;

        EXPECT_NO_THROW(
            coot::reduce::add_amino_hydrogens_gemmi(H_at_name_1, H_at_name_2, first_neighb,
                                        second_neighb_vec, third_neighb_map, bl_amino, *residue_p, st, chain_p)
        );

        int atoms_after = static_cast<int>(residue_p->atoms.size());
        EXPECT_EQ(atoms_after, atoms_before);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}