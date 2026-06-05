#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, dictionary_residue_restraints_t_init) {
    // Load PDB into gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());

    gemmi::Model &model = st.models[0];

    // Case 1: valid residue (first residue of first chain)
    {
        const gemmi::Residue *residue_p = nullptr;
        if (!model.chains.empty()) {
            gemmi::Chain &chain = model.chains[0];
            if (!chain.residues.empty()) {
                residue_p = &chain.residues[0];
            }
        }
        ASSERT_NE(residue_p, nullptr);

        std::string res_name = residue_p->name;
        int seqnum = residue_p->seqid.num.value;

        EXPECT_EQ(res_name, "MET");
        EXPECT_EQ(seqnum, 1);

        // Create an empty dictionary_residue_restraints_t and check state before init
        coot::dictionary_residue_restraints_t drr;

        std::size_t bonds_before = drr.bond_restraint.size();
        std::size_t angles_before = drr.angle_restraint.size();
        std::size_t atoms_before = drr.atom_info.size();

        EXPECT_EQ(bonds_before, 0u);
        EXPECT_EQ(angles_before, 0u);
        EXPECT_EQ(atoms_before, 0u);

        // Call init
        drr.init_gemmi(residue_p);

        // Record state AFTER init
        std::size_t bonds_after = drr.bond_restraint.size();
        std::size_t angles_after = drr.angle_restraint.size();
        std::size_t atoms_after = drr.atom_info.size();

        EXPECT_EQ(bonds_after, 7u);
        EXPECT_EQ(angles_after, 7u);
        EXPECT_EQ(atoms_after, 8u);

        EXPECT_EQ(drr.residue_info.comp_id, "MET");

        // First bond
        ASSERT_FALSE(drr.bond_restraint.empty());
        EXPECT_EQ(drr.bond_restraint[0].atom_id_1(), " CA ");
        EXPECT_EQ(drr.bond_restraint[0].atom_id_2(), " N  ");
        EXPECT_NEAR(drr.bond_restraint[0].dist_, 1.48864, 1e-4);

        // Second bond
        ASSERT_GT(drr.bond_restraint.size(), 1u);
        EXPECT_EQ(drr.bond_restraint[1].atom_id_1(), " CA ");
        EXPECT_EQ(drr.bond_restraint[1].atom_id_2(), " C  ");
        EXPECT_NEAR(drr.bond_restraint[1].dist_, 1.52597, 1e-4);

        // First angle
        ASSERT_FALSE(drr.angle_restraint.empty());
        EXPECT_EQ(drr.angle_restraint[0].atom_id_1(), " N  ");
        EXPECT_EQ(drr.angle_restraint[0].atom_id_2(), " CA ");
        EXPECT_EQ(drr.angle_restraint[0].atom_id_3(), " C  ");
        EXPECT_NEAR(drr.angle_restraint[0].angle_, 109.545, 1e-3);
    }

    // Case 2: null residue — verifies the guarded path (nothing should change)
    {
        coot::dictionary_residue_restraints_t drr2;

        std::size_t bonds_before = drr2.bond_restraint.size();
        std::size_t angles_before = drr2.angle_restraint.size();
        std::size_t atoms_before = drr2.atom_info.size();

        EXPECT_EQ(bonds_before, 0u);
        EXPECT_EQ(angles_before, 0u);
        EXPECT_EQ(atoms_before, 0u);

        drr2.init_gemmi(nullptr);

        std::size_t bonds_after = drr2.bond_restraint.size();
        std::size_t angles_after = drr2.angle_restraint.size();
        std::size_t atoms_after = drr2.atom_info.size();

        EXPECT_EQ(bonds_after, 0u);
        EXPECT_EQ(angles_after, 0u);
        EXPECT_EQ(atoms_after, 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}