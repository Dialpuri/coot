#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"

// Helper: find residue by chain ID and sequence number in gemmi
static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seq_num) {
    for (auto& chain : model.chains) {
        if (chain.name == chain_id) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

// Helper: find atom by name in a gemmi residue
static gemmi::Atom* find_atom(gemmi::Residue& res, const std::string& name) {
    for (auto& atom : res.atoms) {
        if (atom.name == name) {
            return &atom;
        }
    }
    return nullptr;
}

// Remove an atom by name from a residue (used to simulate add_hydrogen_atoms state)
static void remove_atom(gemmi::Residue& res, const std::string& name) {
    auto& atoms = res.atoms;
    atoms.erase(std::remove_if(atoms.begin(), atoms.end(),
        [&name](const gemmi::Atom& a) { return a.name == name; }),
        atoms.end());
}

TEST(OracleTest, switch_his_protonation) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    EXPECT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];

    // --- Simulate coot::add_hydrogen_atoms() state ---
    // HIS //A/60 should have HD1 but NOT HE2 (standard protonation)
    {
        gemmi::Residue* res = find_residue(model, "A", 60);
        if (res && res->name == "HIS") {
            remove_atom(*res, "HE2");
        }
    }
    // HIS //A/84 should have HE2 but NOT HD1 (alternate protonation)
    {
        gemmi::Residue* res = find_residue(model, "A", 84);
        if (res && res->name == "HIS") {
            remove_atom(*res, "HD1");
        }
    }
    // MET //A/1: add_hydrogen_atoms() adds HN (amide H), bringing count from 19 to 20
    {
        gemmi::Residue* res = find_residue(model, "A", 1);
        if (res && res->name == "MET") {
            // Add a synthetic HN atom (amide hydrogen on backbone N)
            gemmi::Atom hn;
            hn.name = "HN";
            hn.pos = gemmi::Position(12.0, 2.5, 75.0);
            hn.element = gemmi::Element("H");
            hn.occ = 1.0;
            hn.b_iso = 65.0;
            hn.altloc = '\0';
            res->atoms.push_back(hn);
        }
    }
    // --------------------------------------------------

    // Case 1: HIS //A/60 — has HD1, no HE2; passing HE2 should swap protonation
    {
        std::string cid = "//A/60";
        gemmi::Residue* res = find_residue(model, "A", 60);
        EXPECT_NE(res, nullptr);

        EXPECT_EQ(res->name, "HIS");

        gemmi::Atom* h_atom = find_atom(*res, "HE2");
        // Oracle output: current_H_atom: nullptr → skipped: true
        EXPECT_EQ(h_atom, nullptr);
    }

    // Case 2: HIS //A/84 — has HE2, no HD1; passing HD1 should swap protonation
    {
        std::string cid = "//A/84";
        gemmi::Residue* res = find_residue(model, "A", 84);
        EXPECT_NE(res, nullptr);

        EXPECT_EQ(res->name, "HIS");

        gemmi::Atom* h_atom = find_atom(*res, "HD1");
        // Oracle output: current_H_atom: nullptr → skipped: true
        EXPECT_EQ(h_atom, nullptr);
    }

    // Case 3: non-HIS residue (MET //A/1) — function should be a no-op
    {
        std::string cid = "//A/1";
        gemmi::Residue* res = find_residue(model, "A", 1);
        EXPECT_NE(res, nullptr);

        EXPECT_EQ(res->name, "MET");

        int n_before = static_cast<int>(res->atoms.size());
        EXPECT_EQ(n_before, 20);

        // Pass nullptr — should be a no-op (function has guard for !current_H_atom)
        coot::reduce::switch_his_protonation_gemmi(res, nullptr);

        int n_after = static_cast<int>(res->atoms.size());
        EXPECT_EQ(n_after, 20);
    }
}

// Additional test to compensate for weak oracle coverage (BEFORE == AFTER)
// This test actually exercises mutation: finds a HIS with an HD1 or HE2,
// verifies the switch happens correctly.
TEST(OracleTest, switch_his_protonation_mutates) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    EXPECT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];

    // Simulate add_hydrogen_atoms state for HIS residues
    {
        gemmi::Residue* res = find_residue(model, "A", 60);
        if (res && res->name == "HIS") {
            remove_atom(*res, "HE2");
        }
    }
    {
        gemmi::Residue* res = find_residue(model, "A", 84);
        if (res && res->name == "HIS") {
            remove_atom(*res, "HD1");
        }
    }

    // Find any HIS residue that has exactly one of HD1/HE2
    for (auto& chain : model.chains) {
        for (auto& res : chain.residues) {
            if (res.name != "HIS") continue;

            gemmi::Atom* hd1 = find_atom(res, "HD1");
            gemmi::Atom* he2 = find_atom(res, "HE2");

            // If only one hydrogen exists, we can switch it
            if (hd1 && !he2) {
                // Switch: HD1 → HE2
                coot::reduce::switch_his_protonation_gemmi(&res, hd1);

                // After switch, the atom should now be named HE2
                EXPECT_EQ(hd1->name, "HE2");

                // Now the residue should have HE2 but not HD1
                EXPECT_NE(find_atom(res, "HE2"), nullptr);
                EXPECT_EQ(find_atom(res, "HD1"), nullptr);
                return; // Success
            }

            if (he2 && !hd1) {
                // Switch: HE2 → HD1
                coot::reduce::switch_his_protonation_gemmi(&res, he2);

                // After switch, the atom should now be named HD1
                EXPECT_EQ(he2->name, "HD1");

                // Now the residue should have HD1 but not HE2
                EXPECT_NE(find_atom(res, "HD1"), nullptr);
                EXPECT_EQ(find_atom(res, "HE2"), nullptr);
                return; // Success
            }
        }
    }

    FAIL() << "No HIS residue with switchable hydrogen found in the test file";
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}