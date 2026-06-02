#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, match_container_t_find_best_match) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model &model = st.models[0];

    // Find chain A
    gemmi::Chain *chainA = nullptr;
    for (gemmi::Chain &ch : model.chains) {
        if (ch.name == "A") { chainA = &ch; break; }
    }
    ASSERT_NE(chainA, nullptr);

    // Find ALA 21 and ALA 31 in chain A
    gemmi::Residue *res_21 = nullptr;
    gemmi::Residue *res_31 = nullptr;
    for (gemmi::Residue &r : chainA->residues) {
        if (r.seqid.num.value == 21) res_21 = &r;
        if (r.seqid.num.value == 31) res_31 = &r;
    }
    ASSERT_NE(res_21, nullptr);
    ASSERT_NE(res_31, nullptr);

    // Collect atoms by name from each residue
    std::map<std::string, gemmi::Atom*> atoms_21;
    std::map<std::string, gemmi::Atom*> atoms_31;

    for (gemmi::Atom &a : res_21->atoms) {
        atoms_21[a.name] = &a;
    }
    for (gemmi::Atom &a : res_31->atoms) {
        atoms_31[a.name] = &a;
    }

    // --- Case 1: Two residues with >2 matching atom pairs ---
    {
        coot::match_container_t match_set;
        int pairs_added = 0;
        for (auto &p : atoms_21) {
            auto it = atoms_31.find(p.first);
            if (it != atoms_31.end()) {
                match_set.add_gemmi(p.second, it->second, model);
                pairs_added++;
            }
        }
        EXPECT_EQ(pairs_added, 5);

        coot::match_container_for_residues_t_gemmi best = match_set.find_best_match_gemmi();

        EXPECT_NE(best.residue_1, nullptr);
        if (best.residue_1) {
            std::string res1_str = best.chain_name_1
                + "/" + std::to_string(best.residue_1->seqid.num.value)
                + " " + best.residue_1->name;
            EXPECT_EQ(res1_str, "A/21 ALA");

            std::string res2_str = best.chain_name_2
                + "/" + std::to_string(best.residue_2->seqid.num.value)
                + " " + best.residue_2->name;
            EXPECT_EQ(res2_str, "A/31 ALA");

            EXPECT_EQ(best.atom_pairs.size(), 5u);
        }
    }

    // --- Case 2: Empty match_container (no matches at all) ---
    {
        coot::match_container_t empty_set;
        coot::match_container_for_residues_t_gemmi empty_best = empty_set.find_best_match_gemmi();
        EXPECT_EQ(empty_best.residue_1, nullptr);
    }

    // --- Case 3: Match with only 2 pairs (< threshold of 3, should return null) ---
    {
        coot::match_container_t small_set;
        gemmi::Atom *n_21  = atoms_21[" N  "];
        gemmi::Atom *ca_21 = atoms_21[" CA "];
        gemmi::Atom *n_31  = atoms_31[" N  "];
        gemmi::Atom *ca_31 = atoms_31[" CA "];
        if (n_21 && ca_21 && n_31 && ca_31) {
            small_set.add_gemmi(n_21, n_31, model);
            small_set.add_gemmi(ca_21, ca_31, model);
        }
        coot::match_container_for_residues_t_gemmi small_best = small_set.find_best_match_gemmi();
        EXPECT_EQ(small_best.residue_1, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}