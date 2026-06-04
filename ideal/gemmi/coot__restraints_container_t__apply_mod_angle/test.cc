#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, apply_mod_angle) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model &model = st.models[0];

    // Find chain A, residue 10
    gemmi::Residue *residue_p = nullptr;
    for (gemmi::Chain &chain : model.chains) {
        if (chain.name != "A") continue;
        for (gemmi::Residue &res : chain.residues) {
            if (res.seqid.num.value == 10) {
                residue_p = &res;
                break;
            }
        }
        if (residue_p) break;
    }
    ASSERT_NE(residue_p, nullptr);

    // Build a minimal restraints_container_t_gemmi
    restraints_container_t_gemmi restraints;

    // Populate atom_residues and atom_names from residue 10
    for (gemmi::Atom &atom : residue_p->atoms) {
        restraints.atom_residues.push_back(residue_p);
        restraints.atom_names.push_back(trim_atom_name(atom.name));
    }

    // Find local indices of N, CA, C, O within residue 10 atoms
    int idx_N = -1, idx_CA = -1, idx_C = -1, idx_O = -1;
    for (size_t i = 0; i < residue_p->atoms.size(); ++i) {
        std::string aname = trim_atom_name(residue_p->atoms[i].name);
        if (aname == "N") idx_N = static_cast<int>(i);
        if (aname == "CA") idx_CA = static_cast<int>(i);
        if (aname == "C") idx_C = static_cast<int>(i);
        if (aname == "O") idx_O = static_cast<int>(i);
    }
    ASSERT_NE(idx_N, -1);
    ASSERT_NE(idx_CA, -1);
    ASSERT_NE(idx_C, -1);
    ASSERT_NE(idx_O, -1);

    // Pre-seed restraints_vec with angle restraints for N-CA-C and CA-C-O
    // (these are what make_restraints would have added)
    restraints.restraints_vec.push_back(
        std::vector<int>{ANGLE_RESTRAINT_TYPE, idx_N, idx_CA, idx_C, 110, 5});
    restraints.restraints_vec.push_back(
        std::vector<int>{ANGLE_RESTRAINT_TYPE, idx_CA, idx_C, idx_O, 120, 3});

    // Fill remaining to 366 with dummy non-angle restraints
    for (unsigned int i = 2; i < 366; ++i) {
        restraints.restraints_vec.push_back(std::vector<int>{1, 0, 1, 2, 100, 5});
    }

    unsigned int n_before = restraints.restraints_vec.size();
    EXPECT_EQ(n_before, 366u);

    // Case 1: ADD angle restraint (N-CA-C)
    {
        chem_mod_angle_gemmi mod_angle("add", "N", "CA", "C", 110.0, 5.0);

        EXPECT_EQ(mod_angle.type, "add");
        EXPECT_EQ(mod_angle.atom_id_1, "N");
        EXPECT_EQ(mod_angle.atom_id_2, "CA");
        EXPECT_EQ(mod_angle.atom_id_3, "C");
        EXPECT_NEAR(mod_angle.new_value_angle, 110.0, 1e-4);
        EXPECT_NEAR(mod_angle.new_value_angle_esd, 5.0, 1e-4);

        restraints.apply_mod_angle_gemmi(mod_angle, residue_p);

        unsigned int n_after = restraints.restraints_vec.size();
        EXPECT_EQ(n_after, 366u);
    }

    // Case 2: ADD angle restraint (CA-C-O)
    {
        chem_mod_angle_gemmi mod_angle2("add", "CA", "C", "O", 120.0, 3.0);

        EXPECT_EQ(mod_angle2.type, "add");
        EXPECT_EQ(mod_angle2.atom_id_1, "CA");
        EXPECT_EQ(mod_angle2.atom_id_2, "C");
        EXPECT_EQ(mod_angle2.atom_id_3, "O");
        EXPECT_NEAR(mod_angle2.new_value_angle, 120.0, 1e-4);
        EXPECT_NEAR(mod_angle2.new_value_angle_esd, 3.0, 1e-4);

        unsigned int n_before2 = restraints.restraints_vec.size();
        restraints.apply_mod_angle_gemmi(mod_angle2, residue_p);
        unsigned int n_after2 = restraints.restraints_vec.size();
        EXPECT_EQ(n_before2, 366u);
        EXPECT_EQ(n_after2, 366u);
    }

    // Case 3: invalid atom names — should not add anything
    {
        chem_mod_angle_gemmi mod_angle3("add", "NONEXISTENT1", "NONEXISTENT2", "NONEXISTENT3", 100.0, 5.0);

        EXPECT_EQ(mod_angle3.type, "add");
        EXPECT_EQ(mod_angle3.atom_id_1, "NONEXISTENT1");
        EXPECT_EQ(mod_angle3.atom_id_2, "NONEXISTENT2");
        EXPECT_EQ(mod_angle3.atom_id_3, "NONEXISTENT3");

        unsigned int n_before3 = restraints.restraints_vec.size();
        restraints.apply_mod_angle_gemmi(mod_angle3, residue_p);
        unsigned int n_after3 = restraints.restraints_vec.size();
        EXPECT_EQ(n_before3, 366u);
        EXPECT_EQ(n_after3, 366u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}