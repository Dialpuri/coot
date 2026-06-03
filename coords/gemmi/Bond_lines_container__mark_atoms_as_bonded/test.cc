#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <map>
#include <string>

static gemmi::Atom* find_atom_by_cid(gemmi::Model& model, const std::string& cid) {
    // cid format: "//CHAIN/SEQNUM/ATMNAME", e.g. "//A/1/N"
    std::string chain_id = cid.substr(2);
    auto pos = chain_id.find('/');
    if (pos == std::string::npos) return nullptr;
    std::string chain_name(1, chain_id[0]);
    std::string rest = chain_id.substr(pos + 1);
    pos = rest.find('/');
    if (pos == std::string::npos) return nullptr;
    int seq_num = std::stoi(rest.substr(0, pos));
    std::string atom_name = rest.substr(pos + 1);

    for (auto& chain : model.chains) {
        if (chain.name == chain_name) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    for (auto& atom : res.atoms) {
                        // atom.name is 4-char padded, e.g. " N  " or " CA "
                        // Compare trimmed
                        std::string trimmed = atom.name;
                        trimmed.erase(0, trimmed.find_first_not_of(" "));
                        trimmed.erase(trimmed.find_last_not_of(" ") + 1);
                        if (trimmed == atom_name) {
                            return &atom;
                        }
                    }
                }
            }
        }
    }
    return nullptr;
}

static int bonded_map_get(const std::map<gemmi::Atom*, int>& bonded_map, gemmi::Atom* atom) {
    auto it = bonded_map.find(atom);
    return (it != bonded_map.end()) ? it->second : 0;
}

static void bonded_map_set(std::map<gemmi::Atom*, int>& bonded_map, gemmi::Atom* atom, int val) {
    bonded_map[atom] = val;
}

TEST(OracleTest, mark_atoms_as_bonded) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Pre-populate bonded_map to simulate MMDB's geometry_init_standard
    // which sets all atoms' UDData to 1 by default.
    std::map<gemmi::Atom*, int> bonded_map;
    for (auto& chain : model.chains) {
        for (auto& res : chain.residues) {
            for (auto& atom : res.atoms) {
                bonded_map[&atom] = 1;
            }
        }
    }

    // ─── Case 1: Two standard atoms (N and CA in residue 1, chain A) ───
    {
        gemmi::Atom *atom1 = find_atom_by_cid(model, "//A/1/N");
        gemmi::Atom *atom2 = find_atom_by_cid(model, "//A/1/CA");

        ASSERT_NE(atom1, nullptr);
        ASSERT_NE(atom2, nullptr);

        EXPECT_EQ(std::string(atom1->element.name()), "N");
        EXPECT_EQ(std::string(atom2->element.name()), "C");

        // Read bonded_map before
        int udd_val_1_before = bonded_map_get(bonded_map, atom1);
        int udd_val_2_before = bonded_map_get(bonded_map, atom2);

        EXPECT_EQ(udd_val_1_before, 1);
        EXPECT_EQ(udd_val_2_before, 1);

        // Call the function
        mark_atoms_as_bonded_gemmi(atom1, atom2, true, bonded_map, false);

        // Read bonded_map after
        int udd_val_1_after = bonded_map_get(bonded_map, atom1);
        int udd_val_2_after = bonded_map_get(bonded_map, atom2);

        EXPECT_EQ(udd_val_1_after, 1);
        EXPECT_EQ(udd_val_2_after, 1);
        EXPECT_EQ(graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND, 1);
    }

    // ─── Case 2: Atom with element "S" (sulfur) — should be skipped ───
    {
        gemmi::Atom *atom1 = find_atom_by_cid(model, "//A/118/SG");
        gemmi::Atom *atom2 = find_atom_by_cid(model, "//A/1/CA");

        ASSERT_NE(atom1, nullptr);
        ASSERT_NE(atom2, nullptr);

        EXPECT_EQ(std::string(atom1->element.name()), "S");
        EXPECT_EQ(std::string(atom2->element.name()), "C");

        int udd_val_1_before = bonded_map_get(bonded_map, atom1);
        int udd_val_2_before = bonded_map_get(bonded_map, atom2);

        EXPECT_EQ(udd_val_1_before, 1);
        EXPECT_EQ(udd_val_2_before, 1);

        mark_atoms_as_bonded_gemmi(atom1, atom2, true, bonded_map, false);

        int udd_val_1_after = bonded_map_get(bonded_map, atom1);
        int udd_val_2_after = bonded_map_get(bonded_map, atom2);

        EXPECT_EQ(udd_val_1_after, 1);
        EXPECT_EQ(udd_val_2_after, 1);
    }

    // ─── Case 3: done_bond_udd_handle=true — should skip ALL writes ───
    {
        gemmi::Atom *atom1 = find_atom_by_cid(model, "//A/2/N");
        gemmi::Atom *atom2 = find_atom_by_cid(model, "//A/2/CA");

        ASSERT_NE(atom1, nullptr);
        ASSERT_NE(atom2, nullptr);

        EXPECT_EQ(std::string(atom1->element.name()), "N");
        EXPECT_EQ(std::string(atom2->element.name()), "C");

        int udd_val_1_before = bonded_map_get(bonded_map, atom1);
        int udd_val_2_before = bonded_map_get(bonded_map, atom2);

        EXPECT_EQ(udd_val_1_before, 1);
        EXPECT_EQ(udd_val_2_before, 1);

        mark_atoms_as_bonded_gemmi(atom1, atom2, true, bonded_map, true);

        int udd_val_1_after = bonded_map_get(bonded_map, atom1);
        int udd_val_2_after = bonded_map_get(bonded_map, atom2);

        EXPECT_EQ(udd_val_1_after, 1);
        EXPECT_EQ(udd_val_2_after, 1);
    }

    // ─── Case 4 (compensation): Set bonded_map to 0 first, then call with standard atoms ───
    // This tests that the function actually writes BONDED_WITH_STANDARD_ATOM_BOND (1)
    {
        gemmi::Atom *atom1 = find_atom_by_cid(model, "//A/3/N");
        gemmi::Atom *atom2 = find_atom_by_cid(model, "//A/3/CA");

        ASSERT_NE(atom1, nullptr);
        ASSERT_NE(atom2, nullptr);

        // Reset bonded_map to 0 to detect actual mutation
        bonded_map_set(bonded_map, atom1, 0);
        bonded_map_set(bonded_map, atom2, 0);

        int udd_val_1_before = bonded_map_get(bonded_map, atom1);
        int udd_val_2_before = bonded_map_get(bonded_map, atom2);

        EXPECT_EQ(udd_val_1_before, 0);
        EXPECT_EQ(udd_val_2_before, 0);

        mark_atoms_as_bonded_gemmi(atom1, atom2, true, bonded_map, false);

        int udd_val_1_after = bonded_map_get(bonded_map, atom1);
        int udd_val_2_after = bonded_map_get(bonded_map, atom2);

        // After calling, atoms should be marked as BONDED_WITH_STANDARD_ATOM_BOND (value 1)
        EXPECT_EQ(udd_val_1_after, 1);
        EXPECT_EQ(udd_val_2_after, 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}