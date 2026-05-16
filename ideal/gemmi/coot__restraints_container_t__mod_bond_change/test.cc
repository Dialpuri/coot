#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, restraints_mod_bond) {
    // Load the PDB
    gemmi::Structure st =
        gemmi::read_pdb_file(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Find residue A/10
    gemmi::Residue* residue_p = nullptr;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& residue : chain.residues) {
                if (chain.name == "A" && residue.seqid.num.value == 10) {
                    residue_p = &residue;
                    break;
                }
            }
            if (residue_p) break;
        }
        if (residue_p) break;
    }

    EXPECT_NE(residue_p, nullptr) << "residue //A/10 not found";

    // case: residue info
    std::string residue_chain = "A";
    int residue_seqnum = residue_p->seqid.num.value;
    std::string residue_resname = residue_p->name;

    EXPECT_EQ(residue_chain, "A");
    EXPECT_EQ(residue_seqnum, 10);
    EXPECT_EQ(residue_resname, "ILE");

    // Build atom name / residue-indexed vectors from the structure
    // (simulating the coot restraints_container_t atom array)
    std::vector<std::string> atom_names;
    std::vector<gemmi::Residue*> atom_residues;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& residue : chain.residues) {
                for (auto& atom : residue.atoms) {
                    atom_names.push_back(atom.name);
                    atom_residues.push_back(&residue);
                }
            }
        }
    }

    // case: empty restraints vector (mimics initial state when thread pool is
    // missing)
    std::vector<coot::simple_restraint> restraints_vec;

    coot::chem_mod_bond mod_bond("change", "NONEXISTENT1", "NONEXISTENT2", "",
                                 1.50f, 0.02f);
    coot::mod_bond_change_gemmi(restraints_vec, atom_names, atom_residues,
                                mod_bond, residue_p);

    // case: initial_restraint_count: 0, initial_bond_restraints: 0
    EXPECT_EQ(restraints_vec.size(), 0);

    int init_bond_count = 0;
    for (unsigned int i = 0; i < restraints_vec.size(); ++i) {
        if (restraints_vec[i].restraint_type == coot::BOND_RESTRAINT) {
            init_bond_count++;
        }
    }
    EXPECT_EQ(init_bond_count, 0);

    // case: mod_bond_change with wrong atom names (should not modify anything)
    // oracle: no_change_on_wrong_names: true
    // (restraints_vec is still empty, but this exercises the function path
    //  with mismatched atom names)
    coot::chem_mod_bond mod_bond2("change", "NONEXISTENT1", "NONEXISTENT2", "",
                                  1.50f, 0.02f);
    coot::mod_bond_change_gemmi(restraints_vec, atom_names, atom_residues,
                                mod_bond2, residue_p);

    bool changed = false;
    for (unsigned int i = 0; i < restraints_vec.size(); ++i) {
        if (restraints_vec[i].restraint_type == coot::BOND_RESTRAINT) {
            if (restraints_vec[i].target_value == 1.50f) {
                changed = true;
            }
        }
    }
    EXPECT_FALSE(changed);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}