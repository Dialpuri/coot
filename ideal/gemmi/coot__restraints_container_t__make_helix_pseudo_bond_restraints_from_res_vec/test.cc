#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find residue in chain by sequence number
static gemmi::Residue* find_residue_by_seq(gemmi::Chain& chain, int seq) {
    for (auto& res : chain.residues) {
        if (res.seqid.num.value == seq)
            return &res;
    }
    return nullptr;
}

TEST(OracleTest, make_helix_pseudo_bond_restraints_from_res_vec) {
    {
        // Case: normal 11-residue segment from chain A (oracle ground truth)
        gemmi::Structure st = gemmi::read_pdb_file(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        std::vector<std::pair<bool, gemmi::CRA>> residues;
        for (int seq = 10; seq <= 20; seq++) {
            for (auto& chain : st.models[0].chains) {
                if (chain.name != "A") continue;
                gemmi::Residue* res = find_residue_by_seq(chain, seq);
                if (res) {
                    gemmi::CRA cra{&chain, res, nullptr};
                    residues.push_back(std::make_pair(false, cra));
                }
            }
        }

        EXPECT_EQ(residues.size(), 11u);

        coot::restraints_container_t restraints;
        restraints.residues_vec = residues;

        size_t restraints_before = restraints.restraints_vec.size();
        EXPECT_EQ(restraints_before, 0u);

        restraints.make_helix_pseudo_bond_restraints_from_res_vec_gemmi();

        size_t restraints_after = restraints.restraints_vec.size();
        EXPECT_EQ(restraints_after, 30u);
        EXPECT_EQ(restraints_after - restraints_before, 30u);
    }

    {
        // Case: small 4-residue segment to exercise fewer branches
        gemmi::Structure st = gemmi::read_pdb_file(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        std::vector<std::pair<bool, gemmi::CRA>> residues;
        for (int seq = 10; seq <= 13; seq++) {
            for (auto& chain : st.models[0].chains) {
                if (chain.name != "A") continue;
                gemmi::Residue* res = find_residue_by_seq(chain, seq);
                if (res) {
                    gemmi::CRA cra{&chain, res, nullptr};
                    residues.push_back(std::make_pair(false, cra));
                }
            }
        }

        EXPECT_EQ(residues.size(), 4u);

        coot::restraints_container_t restraints;
        restraints.residues_vec = residues;

        size_t restraints_before = restraints.restraints_vec.size();
        EXPECT_EQ(restraints_before, 0u);

        restraints.make_helix_pseudo_bond_restraints_from_res_vec_gemmi();

        size_t restraints_after = restraints.restraints_vec.size();
        EXPECT_EQ(restraints_after, 2u);
        EXPECT_EQ(restraints_after - restraints_before, 2u);
    }

    {
        // Case: empty residue vector exercises early return
        gemmi::Structure st = gemmi::read_pdb_file(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        std::vector<std::pair<bool, gemmi::CRA>> residues;

        coot::restraints_container_t restraints;
        restraints.residues_vec = residues;

        size_t restraints_before = restraints.restraints_vec.size();
        EXPECT_EQ(restraints_before, 0u);

        restraints.make_helix_pseudo_bond_restraints_from_res_vec_gemmi();

        size_t restraints_after = restraints.restraints_vec.size();
        EXPECT_EQ(restraints_after, 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}