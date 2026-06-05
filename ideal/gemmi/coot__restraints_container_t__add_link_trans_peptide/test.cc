#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a residue by chain name and residue sequence number in a gemmi Model
static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seq_num) {
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

// Helper: find chain by name
static gemmi::Chain* find_chain(gemmi::Model& model, const std::string& chain_id) {
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            return &chain;
        }
    }
    return nullptr;
}

// Helper: build cras vector from all atoms in a model
static std::vector<gemmi::CRA> build_cras(gemmi::Model& model) {
    std::vector<gemmi::CRA> cras;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            for (gemmi::Atom& atom : res.atoms) {
                cras.push_back(gemmi::CRA{&chain, &res, &atom});
            }
        }
    }
    return cras;
}

TEST(OracleTest, add_link_trans_peptide) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // Get two adjacent residues with CA, C, N atoms (MET 1 and GLU 2, chain A)
    gemmi::Residue* res1 = find_residue(model, "A", 1);
    gemmi::Residue* res2 = find_residue(model, "A", 2);

    ASSERT_NE(res1, nullptr);
    ASSERT_NE(res2, nullptr);

    gemmi::Chain* chain_A = find_chain(model, "A");
    ASSERT_NE(chain_A, nullptr);

    // Build cras for atom index lookup
    std::vector<gemmi::CRA> cras = build_cras(model);

    restraints_container_t restraints;

    // Case 1: valid adjacent residues (MET 1 -> GLU 2), add_even_if_cis=false
    {
        size_t before = restraints.restraints_vec.size();
        int n_added = restraints.add_link_trans_peptide_gemmi(
            chain_A, res1, chain_A, res2, false, false, false, cras);
        size_t after = restraints.restraints_vec.size();

        EXPECT_EQ(n_added, 1);
        EXPECT_EQ(before, static_cast<size_t>(0));
        EXPECT_EQ(after, static_cast<size_t>(1));
    }

    // Case 2: same pair with add_even_if_cis=true
    {
        size_t before = restraints.restraints_vec.size();
        int n_added = restraints.add_link_trans_peptide_gemmi(
            chain_A, res1, chain_A, res2, false, false, true, cras);
        size_t after = restraints.restraints_vec.size();

        EXPECT_EQ(n_added, 1);
        EXPECT_EQ(before, static_cast<size_t>(1));
        EXPECT_EQ(after, static_cast<size_t>(2));
    }

    // Case 3: non-adjacent residues (should not add restraint because N-C distance > 2.0)
    {
        gemmi::Residue* res_far1 = find_residue(model, "A", 1);
        gemmi::Residue* res_far2 = find_residue(model, "A", 50);

        ASSERT_NE(res_far1, nullptr);
        ASSERT_NE(res_far2, nullptr);

        size_t before = restraints.restraints_vec.size();
        int n_added = restraints.add_link_trans_peptide_gemmi(
            chain_A, res_far1, chain_A, res_far2, false, false, false, cras);
        size_t after = restraints.restraints_vec.size();

        EXPECT_EQ(n_added, 0);
        EXPECT_EQ(before, static_cast<size_t>(2));
        EXPECT_EQ(after, static_cast<size_t>(2));
    }

    // Case 4: complementary - test with is_fixed_first=true (should still add restraint
    // because distance is valid for adjacent residues)
    {
        // Use another adjacent pair to have a fresh test - residues 49 and 50
        gemmi::Residue* res_49 = find_residue(model, "A", 49);
        gemmi::Residue* res_50 = find_residue(model, "A", 50);

        ASSERT_NE(res_49, nullptr);
        ASSERT_NE(res_50, nullptr);

        size_t before = restraints.restraints_vec.size();
        int n_added = restraints.add_link_trans_peptide_gemmi(
            chain_A, res_49, chain_A, res_50, true, false, false, cras);
        size_t after = restraints.restraints_vec.size();

        EXPECT_EQ(n_added, 1);
        EXPECT_EQ(after, before + 1);
    }

    // Case 5: complementary - test reversed order (should not add because residue1 N and residue2 C
    // won't be adjacent in correct peptide direction)
    {
        size_t before = restraints.restraints_vec.size();
        int n_added = restraints.add_link_trans_peptide_gemmi(
            chain_A, res2, chain_A, res1, false, false, false, cras);
        size_t after = restraints.restraints_vec.size();

        // Reversed residues should not form a trans peptide link
        EXPECT_EQ(n_added, 0);
        EXPECT_EQ(after, before);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}