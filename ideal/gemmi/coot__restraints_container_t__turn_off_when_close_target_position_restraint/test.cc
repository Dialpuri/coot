#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTestGemmi, turn_off_when_close_target_position_restraint) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model *model = &st.models[0];
    ASSERT_FALSE(model->chains.empty());

    gemmi::Chain *chain = &model->chains[0];
    ASSERT_FALSE(chain->residues.empty());

    gemmi::Residue *residue = &chain->residues[0];
    ASSERT_FALSE(residue->atoms.empty());

    // Build restraints container
    restraints_container_t_gemmi restraints;

    // Collect all atoms from the residue into restraints.atom
    for (const auto &atom : residue->atoms) {
        restraints.atom.push_back(&atom);
    }

    gemmi::Atom *atom = &residue->atoms[0];
    ASSERT_NE(atom, nullptr);

    // Find the atom index in restraints.atom vector
    int atom_index = -1;
    for (int i = 0; i < static_cast<int>(restraints.atom.size()); i++) {
        if (restraints.atom[i] == atom) {
            atom_index = i;
            break;
        }
    }
    ASSERT_GE(atom_index, 0);

    // CASE 1: Target position very close to atom (within 0.6 Å), is_closed = true
    {
        gemmi::Vec3 close_target(atom->pos.x + 0.1, atom->pos.y + 0.0, atom->pos.z + 0.0);

        simple_restraint_gemmi r;
        r.restraint_type = TARGET_POS_RESTRAINT;
        r.atom_index_1 = atom_index;
        r.atom_pull_target_pos = close_target;
        r.is_closed = true;
        restraints.restraints_vec.push_back(r);

        bool result = restraints.turn_off_when_close_target_position_restraint_gemmi();

        EXPECT_TRUE(result);
    }

    // CASE 2: Target position far from atom (> 0.6 Å), is_closed = true
    {
        gemmi::Vec3 far_target(atom->pos.x + 10.0, atom->pos.y + 10.0, atom->pos.z + 10.0);

        simple_restraint_gemmi r;
        r.restraint_type = TARGET_POS_RESTRAINT;
        r.atom_index_1 = atom_index;
        r.atom_pull_target_pos = far_target;
        r.is_closed = true;
        restraints.restraints_vec.push_back(r);

        bool result2 = restraints.turn_off_when_close_target_position_restraint_gemmi();

        EXPECT_TRUE(result2);
    }

    // CASE 3: is_closed = false — should not be actioned even if close
    {
        gemmi::Vec3 close_target3(atom->pos.x + 0.05, atom->pos.y + 0.0, atom->pos.z + 0.0);

        simple_restraint_gemmi r;
        r.restraint_type = TARGET_POS_RESTRAINT;
        r.atom_index_1 = atom_index;
        r.atom_pull_target_pos = close_target3;
        r.is_closed = false;  // Leave is_closed = false (default)
        restraints.restraints_vec.push_back(r);

        bool result3 = restraints.turn_off_when_close_target_position_restraint_gemmi();

        EXPECT_TRUE(result3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}