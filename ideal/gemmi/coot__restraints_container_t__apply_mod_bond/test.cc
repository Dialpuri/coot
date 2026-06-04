#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

namespace coot {

static gemmi::Residue* find_residue(gemmi::Model& model,
                                     const std::string& chain_name,
                                     int seq_num) {
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_name) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seq_num)
                    return &res;
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, apply_mod_bond) {
    // --- Setup: load PDB, select residue ---
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);

    // Get residue using gemmi: chain A, residue 1
    gemmi::Residue* residue_p = find_residue(st.models[0], "A", 1);
    EXPECT_NE(residue_p, nullptr);

    // Verify residue identity matches oracle input
    EXPECT_EQ(residue_p->name, "MET");

    // --- Build restraints vector ---
    std::vector<simple_restraint> restraints_vec;

    int n_restraints_before = static_cast<int>(restraints_vec.size());
    EXPECT_EQ(n_restraints_before, 0);

    // --- Case 1: ADD a bond between CA and CB in this residue ---
    {
        coot::chem_mod_bond mod_bond("add",
                                     "CA",
                                     "CB",
                                     "",
                                     1.54,
                                     0.02);

        EXPECT_EQ(mod_bond.atom_id_1, "CA");
        EXPECT_EQ(mod_bond.atom_id_2, "CB");
        EXPECT_NEAR(mod_bond.new_value_dist, 1.54, 1e-4);
        EXPECT_NEAR(mod_bond.new_value_dist_esd, 0.02, 1e-4);

        apply_mod_bond_gemmi(restraints_vec, mod_bond, *residue_p);

        int n_restraints_after = static_cast<int>(restraints_vec.size());
        EXPECT_EQ(n_restraints_after - n_restraints_before, 1);
    }

    // --- Case 2: DELETE a bond where atom2 does not exist (no-op expected) ---
    {
        coot::chem_mod_bond mod_bond("delete",
                                     "CA",
                                     "NONEXISTENT",
                                     "",
                                     1.54,
                                     0.02);

        EXPECT_EQ(mod_bond.atom_id_1, "CA");
        EXPECT_EQ(mod_bond.atom_id_2, "NONEXISTENT");

        int n_before_del = static_cast<int>(restraints_vec.size());
        apply_mod_bond_gemmi(restraints_vec, mod_bond, *residue_p);
        int n_after_del = static_cast<int>(restraints_vec.size());

        EXPECT_EQ(n_before_del - n_after_del, 0);
    }
}

} // namespace coot

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}