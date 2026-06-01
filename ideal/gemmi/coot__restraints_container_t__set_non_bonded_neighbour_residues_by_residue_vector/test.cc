#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static std::string get_chain_id(const gemmi::Model &model, const gemmi::Residue *res) {
    for (const gemmi::Chain &chain : model.chains) {
        for (const gemmi::Residue &r : chain.residues) {
            if (&r == res) return chain.name;
        }
    }
    return "";
}

TEST(OracleTest_gemmi, set_non_bonded_neighbour_residues_by_residue_vector) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model &model = st.models[0];
    gemmi::Chain &chain_a = model.chains[0];

    gemmi::Residue *res_a10 = &chain_a.residues[9];   // seqnum 10
    gemmi::Residue *res_a11 = &chain_a.residues[10];  // seqnum 11
    gemmi::Residue *res_a5  = &chain_a.residues[4];   // seqnum 5
    gemmi::Residue *res_a15 = &chain_a.residues[14];  // seqnum 15
    gemmi::Residue *res_a20 = &chain_a.residues[19];  // seqnum 20

    ASSERT_NE(res_a10, nullptr);
    ASSERT_NE(res_a11, nullptr);
    ASSERT_NE(res_a5, nullptr);
    ASSERT_NE(res_a15, nullptr);
    ASSERT_NE(res_a20, nullptr);

    // Build residues_vec: moving residues (first=false)
    std::vector<std::pair<bool, gemmi::Residue *>> residues_vec;
    residues_vec.push_back({false, res_a10});  // moving: A/10 ALA
    residues_vec.push_back({false, res_a11});  // moving: A/11 SER

    EXPECT_EQ(residues_vec.size(), 2u);

    coot::protein_geometry pg;
    pg.init_standard();

    // Case 1: neighbours NOT in residues_vec, no bonded flankers
    {
        std::map<gemmi::Residue *, std::set<gemmi::Residue *>> neighbour_set;
        neighbour_set[res_a10].insert(res_a5);   // A/5 is a neighbour of A/10
        neighbour_set[res_a11].insert(res_a15);  // A/15 is a neighbour of A/11

        EXPECT_EQ(neighbour_set.size(), 2u);

        std::vector<coot::bonded_pair_gemmi> bonded_flanking_pairs;
        EXPECT_EQ(bonded_flanking_pairs.size(), 0u);

        // The original test expected before_size=9u (from the existing restraints_container_t state),
        // but in the gemmi port the function is a free function that returns a new vector,
        // so the "before" concept doesn't apply. We verify the result matches expected.

        std::vector<gemmi::Residue *> result =
            coot::set_non_bonded_neighbour_residues_by_residue_vector_gemmi(
                neighbour_set, residues_vec, bonded_flanking_pairs, pg);

        size_t after_size = result.size();
        EXPECT_EQ(after_size, 2u);

        // Check the neighbour residues
        gemmi::Residue *r0 = result[0];
        EXPECT_EQ(get_chain_id(model, r0), "A");
        EXPECT_EQ(r0->seqid.num.value, 5);
        EXPECT_EQ(std::string(r0->name), "GLN");

        gemmi::Residue *r1 = result[1];
        EXPECT_EQ(get_chain_id(model, r1), "A");
        EXPECT_EQ(r1->seqid.num.value, 15);
        EXPECT_EQ(std::string(r1->name), "TYR");
    }

    // Case 2: neighbours already in residues_vec (should be excluded)
    {
        std::map<gemmi::Residue *, std::set<gemmi::Residue *>> neighbour_set2;
        neighbour_set2[res_a10].insert(res_a11); // A/11 is in residues_vec (moving), should be excluded
        neighbour_set2[res_a10].insert(res_a10); // A/10 is in residues_vec (moving), should be excluded

        EXPECT_EQ(neighbour_set2.size(), 1u);

        std::vector<coot::bonded_pair_gemmi> bonded_flanking_pairs;

        std::vector<gemmi::Residue *> result =
            coot::set_non_bonded_neighbour_residues_by_residue_vector_gemmi(
                neighbour_set2, residues_vec, bonded_flanking_pairs, pg);

        size_t after_size2 = result.size();
        EXPECT_EQ(after_size2, 0u);
    }

    // Case 3: some neighbours are bonded flankers (should be excluded)
    {
        std::map<gemmi::Residue *, std::set<gemmi::Residue *>> neighbour_set3;
        neighbour_set3[res_a10].insert(res_a20); // A/20 is a neighbour
        neighbour_set3[res_a10].insert(res_a5);  // A/5 is a neighbour

        EXPECT_EQ(neighbour_set3.size(), 1u);

        std::vector<coot::bonded_pair_gemmi> bonded_pairs_with_flankers;
        // A/20 is in bonded_flanking_pairs as res_1, so should be excluded
        // A/5 is not in bonded_pairs, so should remain
        bonded_pairs_with_flankers.push_back(coot::bonded_pair_gemmi(res_a20, res_a10, false, false, "PEPTIDE"));

        EXPECT_EQ(bonded_pairs_with_flankers.size(), 1u);

        std::vector<gemmi::Residue *> result =
            coot::set_non_bonded_neighbour_residues_by_residue_vector_gemmi(
                neighbour_set3, residues_vec, bonded_pairs_with_flankers, pg);

        size_t after_size3 = result.size();
        EXPECT_EQ(after_size3, 1u);

        // Check the remaining neighbour residue
        gemmi::Residue *r0 = result[0];
        EXPECT_EQ(get_chain_id(model, r0), "A");
        EXPECT_EQ(r0->seqid.num.value, 5);
        EXPECT_EQ(std::string(r0->name), "GLN");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}