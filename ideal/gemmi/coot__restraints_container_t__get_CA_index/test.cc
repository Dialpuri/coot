#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper to find a residue by chain ID and sequence number in a gemmi model
static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seqnum) {
    for (auto& chain : model.chains) {
        if (chain.name == chain_id) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == seqnum) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, restraints_container_t_get_CA_index_gemmi) {
    // Setup: load PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_GE(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];

    clipper::Xmap<float> dummy_xmap;

    // Get residues via gemmi traversal
    // Case 1: residue 10 ILE chain A — has CA
    gemmi::Residue* residue_with_ca = find_residue(model, "A", 10);
    EXPECT_NE(residue_with_ca, nullptr);

    // Case 2: residue 11 GLY chain A — has CA
    gemmi::Residue* residue_gly = find_residue(model, "A", 11);

    std::vector<std::pair<bool, gemmi::Residue*>> residues;
    residues.push_back(std::make_pair(false, residue_with_ca));
    if (residue_gly) {
        residues.push_back(std::make_pair(false, residue_gly));
    }

    // No make_restraints call needed — the gemmi version doesn't need it for get_CA_index
    coot::restraints_container_t_gemmi restraints(residues, "standard", &st, &dummy_xmap);

    // Case 1: residue 10 ILE in chain A — has CA
    {
        int ca_idx = restraints.get_CA_index_gemmi(residue_with_ca);
        EXPECT_EQ(residue_with_ca->name, "ILE");
        EXPECT_EQ(residue_with_ca->seqid.num.value, 10);
        EXPECT_EQ(ca_idx, 1);
    }

    // Case 2: residue 11 GLY in chain A — has CA
    if (residue_gly) {
        int ca_idx = restraints.get_CA_index_gemmi(residue_gly);
        EXPECT_EQ(residue_gly->name, "GLY");
        EXPECT_EQ(residue_gly->seqid.num.value, 11);
        EXPECT_EQ(ca_idx, 1);
    }

    // Case 3: residue NOT in the restraints container — oracle returned 1
    gemmi::Residue* residue_not_in_container = find_residue(model, "A", 20);
    if (residue_not_in_container) {
        int ca_idx = restraints.get_CA_index_gemmi(residue_not_in_container);
        EXPECT_EQ(residue_not_in_container->name, "LYS");
        EXPECT_EQ(residue_not_in_container->seqid.num.value, 20);
        EXPECT_EQ(ca_idx, 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}