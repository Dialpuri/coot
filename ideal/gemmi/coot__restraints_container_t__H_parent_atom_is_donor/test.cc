#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, H_parent_atom_is_donor) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model* model = &st.models[0];
    ASSERT_NE(model, nullptr);

    gemmi::Chain* chain = &model->chains[0];
    ASSERT_NE(chain, nullptr);

    gemmi::Residue* res = &chain->residues[0];
    ASSERT_NE(res, nullptr);

    // Construct restraints_container_t_gemmi minimally
    coot::restraints_container_t_gemmi restraints;

    // Get atoms from the residue
    gemmi::Atom* n_atom = nullptr;
    gemmi::Atom* c_atom = nullptr;
    gemmi::Atom* ca_atom = nullptr;

    for (auto& a : res->atoms) {
        std::string name = a.name;
        if (name == "N" || name == " N  ") n_atom = &a;
        if (name == "C" || name == " C  ") c_atom = &a;
        if (name == "CA" || name == " CA ") ca_atom = &a;
    }

    ASSERT_NE(n_atom, nullptr);
    ASSERT_NE(c_atom, nullptr);
    ASSERT_NE(ca_atom, nullptr);

    // Case 1: N marked as HB_DONOR — should return true
    {
        restraints.H_atom_parent_energy_type_atom_map[n_atom] = coot::hb_t::HB_DONOR;
        bool result = restraints.H_parent_atom_is_donor_gemmi(n_atom);
        EXPECT_TRUE(result);
    }

    // Case 2: C marked as HB_BOTH — should return true
    {
        restraints.H_atom_parent_energy_type_atom_map[c_atom] = coot::hb_t::HB_BOTH;
        bool result = restraints.H_parent_atom_is_donor_gemmi(c_atom);
        EXPECT_TRUE(result);
    }

    // Case 3: CA not in map — should return false (not found)
    {
        bool result = restraints.H_parent_atom_is_donor_gemmi(ca_atom);
        EXPECT_FALSE(result);
    }

    // Case 4: N marked as HB_ACCEPTOR — should return false (acceptor is not donor)
    {
        restraints.H_atom_parent_energy_type_atom_map[n_atom] = coot::hb_t::HB_ACCEPTOR;
        bool result = restraints.H_parent_atom_is_donor_gemmi(n_atom);
        EXPECT_FALSE(result);
    }

    // Case 5: nullptr — should return false
    {
        bool result = restraints.H_parent_atom_is_donor_gemmi(nullptr);
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}