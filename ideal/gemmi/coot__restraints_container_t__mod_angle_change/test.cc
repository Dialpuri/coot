#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, mod_angle_change) {
    // Setup: load PDB, get residue
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model &model = st.models[0];

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

    std::string resname = residue_p->name;
    int seqnum = residue_p->seqid.num.value;

    EXPECT_EQ(resname, "ILE");
    EXPECT_EQ(seqnum, 10);

    // Build a minimal restraints_container_t_gemmi
    restraints_container_t_gemmi restraints;

    // Populate atom_residues and atom_names from residue 10
    for (gemmi::Atom &atom : residue_p->atoms) {
        restraints.atom_residues.push_back(residue_p);
        restraints.atom_names.push_back(trim_atom_name(atom.name));
    }

    // Find local indices of N, CA, C within residue 10
    int idx_N = -1, idx_CA = -1, idx_C = -1;
    for (size_t i = 0; i < residue_p->atoms.size(); ++i) {
        std::string aname = trim_atom_name(residue_p->atoms[i].name);
        if (aname == "N")  idx_N  = static_cast<int>(i);
        if (aname == "CA") idx_CA = static_cast<int>(i);
        if (aname == "C")  idx_C  = static_cast<int>(i);
    }
    ASSERT_NE(idx_N, -1);
    ASSERT_NE(idx_CA, -1);
    ASSERT_NE(idx_C, -1);

    // Pre-seed the N-CA-C angle restraint with the values that
    // make_restraints (via protein_geometry) would produce
    {
        simple_restraint_gemmi rest;
        rest.restraint_type = ANGLE_RESTRAINT;
        rest.atom_index_1   = idx_N;
        rest.atom_index_2   = idx_CA;
        rest.atom_index_3   = idx_C;
        rest.target_value   = 108.738;  // oracle value from original test
        rest.sigma          = 1.5;       // oracle value from original test
        restraints.restraints_vec.push_back(rest);
    }

    // Find the existing ANGLE_RESTRAINT for N-CA-C in this residue
    int found_index = -1;
    double old_target = 0;
    double old_sigma = 0;

    for (unsigned int i = 0; i < restraints.restraints_vec.size(); i++) {
        simple_restraint_gemmi &rest = restraints.restraints_vec[i];
        if (rest.restraint_type == ANGLE_RESTRAINT) {
            if (restraints.atom_residues[static_cast<size_t>(rest.atom_index_1)] == residue_p) {
                if (restraints.atom_residues[static_cast<size_t>(rest.atom_index_2)] == residue_p) {
                    std::string n1 = restraints.atom_names[static_cast<size_t>(rest.atom_index_1)];
                    std::string n2 = restraints.atom_names[static_cast<size_t>(rest.atom_index_2)];
                    std::string n3 = restraints.atom_names[static_cast<size_t>(rest.atom_index_3)];
                    if (n1 == "N" && n2 == "CA" && n3 == "C") {
                        found_index = static_cast<int>(i);
                        old_target = rest.target_value;
                        old_sigma = rest.sigma;
                        break;
                    }
                }
            }
        }
    }

    ASSERT_GE(found_index, 0);

    // Case 1: Oracle output values for old_target and old_sigma
    EXPECT_NEAR(old_target, 108.738, 1e-4);
    EXPECT_NEAR(old_sigma, 1.5, 1e-4);

    // Create the chem_mod_angle with function="change" (uses trimmed names)
    chem_mod_angle_gemmi mod_angle("change", "N", "CA", "C", 120.0, 2.0);

    // Call the gemmi ported function
    restraints.mod_angle_change_gemmi(mod_angle, residue_p);

    double new_target = restraints.restraints_vec[found_index].target_value;
    double new_sigma  = restraints.restraints_vec[found_index].sigma;

    EXPECT_NEAR(new_target, 120.0, 1e-4);
    EXPECT_NEAR(new_sigma, 2.0, 1e-4);

    // Case 2: non-matching angle (should not change anything)
    chem_mod_angle_gemmi mod_angle_bad("change", "OX1", "OX2", "OX3", 100.0, 1.0);

    double target_before_bad = restraints.restraints_vec[found_index].target_value;
    restraints.mod_angle_change_gemmi(mod_angle_bad, residue_p);
    double target_after_bad = restraints.restraints_vec[found_index].target_value;

    EXPECT_TRUE(target_before_bad == target_after_bad);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}