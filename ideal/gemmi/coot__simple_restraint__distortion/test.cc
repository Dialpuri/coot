#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static std::string trim_atom(const std::string &s) {
    return std::string(s).erase(s.find_last_not_of(" \t\r\n") + 1).erase(0, s.find_first_not_of(" \t\r\n"));
}

TEST(OracleTest, simple_restraint_distortion) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model *model = &st.models[0];
    ASSERT_NE(model, nullptr);

    gemmi::Chain *chain = &model->chains[0];
    gemmi::Residue *res = &chain->residues[0];

    std::vector<gemmi::Atom*> atoms;
    for (auto &a : res->atoms) {
        atoms.push_back(&a);
    }

    int natoms = static_cast<int>(atoms.size());
    ASSERT_GE(natoms, 4);

    int n_idx = -1, ca_idx = -1, c_idx = -1, cb_idx = -1;
    for (int i = 0; i < natoms; i++) {
        std::string aname = trim_atom(atoms[i]->name);
        if (aname == "N")  n_idx = i;
        if (aname == "CA") ca_idx = i;
        if (aname == "C")  c_idx = i;
        if (aname == "CB") cb_idx = i;
    }
    ASSERT_NE(n_idx, -1);
    ASSERT_NE(ca_idx, -1);
    ASSERT_NE(c_idx, -1);
    ASSERT_NE(cb_idx, -1);

    // Case 1: BOND_RESTRAINT between N and CA
    {
        coot::simple_restraint rest;
        rest.restraint_type = BOND_RESTRAINT;
        rest.atom_index_1 = n_idx;
        rest.atom_index_2 = ca_idx;
        rest.target_value = 1.45;
        rest.sigma = 0.02;

        std::pair<double, double> result = rest.distortion_gemmi(atoms, 0.0);

        EXPECT_NEAR(result.first, 3.73314, 1e-4);
        EXPECT_NEAR(result.second, 0.0386427, 1e-4);
    }

    // Case 2: ANGLE_RESTRAINT N-CA-C
    {
        coot::simple_restraint rest;
        rest.restraint_type = ANGLE_RESTRAINT;
        rest.atom_index_1 = n_idx;
        rest.atom_index_2 = ca_idx;
        rest.atom_index_3 = c_idx;
        rest.target_value = 110.0;
        rest.sigma = 2.0;

        std::pair<double, double> result = rest.distortion_gemmi(atoms, 0.0);

        EXPECT_NEAR(result.first, 0.0517027, 1e-4);
        EXPECT_NEAR(result.second, -0.454764, 1e-4);
    }

    // Case 3: CHIRAL_VOLUME_RESTRAINT on MET CA (centre=N, C, CB)
    {
        coot::simple_restraint rest;
        rest.restraint_type = CHIRAL_VOLUME_RESTRAINT;
        rest.atom_index_centre = ca_idx;
        rest.atom_index_1 = n_idx;
        rest.atom_index_2 = c_idx;
        rest.atom_index_3 = cb_idx;
        rest.target_chiral_volume = -2.5;
        rest.sigma = 0.2;

        std::pair<double, double> result = rest.distortion_gemmi(atoms, 0.0);

        EXPECT_NEAR(result.first, 653.74, 0.01);
        EXPECT_NEAR(result.second, 5.11367, 1e-4);
    }

    // Case 4: NON_BONDED_CONTACT_RESTRAINT between N and CB
    {
        coot::simple_restraint rest;
        rest.restraint_type = NON_BONDED_CONTACT_RESTRAINT;
        rest.atom_index_1 = n_idx;
        rest.atom_index_2 = cb_idx;
        rest.target_value = 3.4;
        rest.sigma = 0.5;

        std::pair<double, double> result = rest.distortion_gemmi(atoms, 0.1);

        EXPECT_NEAR(result.first, 13.9406, 1e-4);
        EXPECT_NEAR(result.second, -0.906361, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}