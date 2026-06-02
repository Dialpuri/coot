#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTestGemmi, make_helix_pseudo_bond_restraints) {
    // Load PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "PDB load failed";

    gemmi::Model* model = &st.models[0];

    // Original test: mol->GetResidue("//A/" + i) returns NULL for all i=1..10
    // so residues ends up empty. We replicate that.
    std::vector<gemmi::Residue*> residues;
    // (empty — matches original MMDB behavior where GetResidue returns NULL)

    coot::restraints_container_t restraints(model);
    restraints.residues_ = residues;

    // Case 1: from_residue_vector = true (default after construction)
    {
        EXPECT_EQ(static_cast<int>(residues.size()), 0);
        int before = restraints.size();
        restraints.make_helix_pseudo_bond_restraints_gemmi();
        int after = restraints.size();
        EXPECT_EQ(before, 0);
        EXPECT_EQ(after, 0);
    }

    // Case 2: from_residue_vector = false, use chain/residue range
    {
        restraints.from_residue_vector = false;
        restraints.chain_id_save = "A";
        restraints.istart_res = 1;
        restraints.iend_res = 30;

        int before = restraints.size();
        restraints.make_helix_pseudo_bond_restraints_gemmi();
        int after = restraints.size();
        EXPECT_EQ(before, 0);
        EXPECT_EQ(after, 52);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}