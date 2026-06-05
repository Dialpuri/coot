#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, deep_copy_this_residue_add_chain) {
    // Setup: read PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];

    // --- Case 1: valid residue, whole_residue=true, attach_to_new_chain=true ---
    {
        // Find chain A, residue 10 (ILE)
        const gemmi::Chain* chain_A = model.find_chain("A");
        ASSERT_NE(chain_A, nullptr);
        const gemmi::Residue* res = nullptr;
        for (const auto& r : chain_A->residues) {
            if (r.seqid.num.value == 10) {
                res = &r;
                break;
            }
        }
        ASSERT_NE(res, nullptr);

        std::string altconf("");
        bool whole_residue_flag = true;
        bool attach_to_new_chain_flag = true;

        EXPECT_STREQ(res->name.c_str(), "ILE");
        EXPECT_EQ(res->seqid.num.value, 10);

        int atoms_before = static_cast<int>(res->atoms.size());

        gemmi::CRA cra{const_cast<gemmi::Chain*>(chain_A),
                       const_cast<gemmi::Residue*>(res),
                       nullptr};

        auto [copy, chain_p] = coot::util::deep_copy_this_residue_add_chain_gemmi(
            cra, altconf, whole_residue_flag, attach_to_new_chain_flag);

        EXPECT_NE(copy, nullptr);
        if (copy) {
            int atoms_after = static_cast<int>(copy->atoms.size());
            EXPECT_STREQ(copy->name.c_str(), "ILE");
            EXPECT_EQ(copy->seqid.num.value, 10);
            EXPECT_EQ(atoms_before, 8);
            EXPECT_EQ(atoms_after, 8);

            if (chain_p) {
                EXPECT_STREQ(chain_p->name.c_str(), "A");
                delete chain_p; // clean up chain which owns the residue
            } else {
                delete copy;
            }
        }
    }

    // --- Case 2: invalid residue (NULL input) — verifies guarded path ---
    {
        // Find residue 9999 — should not exist
        const gemmi::Chain* chain_A = model.find_chain("A");
        const gemmi::Residue* res = nullptr;
        for (const auto& r : chain_A->residues) {
            if (r.seqid.num.value == 9999) {
                res = &r;
                break;
            }
        }
        EXPECT_EQ(res, nullptr);

        gemmi::CRA cra{const_cast<gemmi::Chain*>(chain_A), nullptr, nullptr};

        auto [copy, chain_p] = coot::util::deep_copy_this_residue_add_chain_gemmi(
            cra, "", true, true);
        EXPECT_EQ(copy, nullptr);
        if (chain_p) {
            delete chain_p;
        }
    }

    // --- Case 3: valid residue, whole_residue=false, specific altconf filter ---
    {
        // Find chain A, residue 50 (ARG)
        const gemmi::Chain* chain_A = model.find_chain("A");
        ASSERT_NE(chain_A, nullptr);
        const gemmi::Residue* res = nullptr;
        for (const auto& r : chain_A->residues) {
            if (r.seqid.num.value == 50) {
                res = &r;
                break;
            }
        }
        ASSERT_NE(res, nullptr);

        std::string altconf("A");
        bool whole_residue_flag = false;
        bool attach_to_new_chain_flag = false;

        EXPECT_STREQ(res->name.c_str(), "ARG");

        int atoms_before = static_cast<int>(res->atoms.size());

        gemmi::CRA cra{const_cast<gemmi::Chain*>(chain_A),
                       const_cast<gemmi::Residue*>(res),
                       nullptr};

        auto [copy, chain_p] = coot::util::deep_copy_this_residue_add_chain_gemmi(
            cra, altconf, whole_residue_flag, attach_to_new_chain_flag);

        EXPECT_NE(copy, nullptr);
        if (copy) {
            int atoms_after = static_cast<int>(copy->atoms.size());
            EXPECT_STREQ(copy->name.c_str(), "ARG");
            EXPECT_EQ(atoms_before, 11);
            EXPECT_EQ(atoms_after, 11);
            delete copy;
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}