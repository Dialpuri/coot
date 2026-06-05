#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, water_coordination_t_transform_atom) {
    // case: identity transform
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];
        gemmi::Residue& res = chain.residues[0];
        gemmi::Atom* atom1 = &res.atoms[0];

        int natoms = res.atoms.size();
        gemmi::Atom* atom2 = nullptr;
        if (natoms > 1) {
            atom2 = &res.atoms[1];
        } else {
            gemmi::Residue* res2 = nullptr;
            if (chain.residues.size() > 1) {
                res2 = &chain.residues[1];
            }
            if (res2 && !res2->atoms.empty()) {
                atom2 = &res2->atoms[0];
            }
        }
        if (!atom2) atom2 = atom1;

        coot::util::water_coordination_t wc;
        coot::util::contact_atoms_info_t cai(atom1, coot::util::contact_atom_t(atom2, atom1));
        wc.atom_contacts.push_back(std::move(cai));

        wc.transform_atom_gemmi(0, 0);

        EXPECT_NEAR(atom2->pos.x, 13.401, 1e-4);
        EXPECT_NEAR(atom2->pos.y, 3.78, 1e-4);
        EXPECT_NEAR(atom2->pos.z, 74.703, 1e-4);
    }

    // case: translation transform
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];
        gemmi::Residue& res = chain.residues[0];
        gemmi::Atom* atom1 = &res.atoms[0];

        int natoms = res.atoms.size();
        gemmi::Atom* atom2 = nullptr;
        if (natoms > 1) {
            atom2 = &res.atoms[1];
        } else {
            gemmi::Residue* res2 = nullptr;
            if (chain.residues.size() > 1) {
                res2 = &chain.residues[1];
            }
            if (res2 && !res2->atoms.empty()) {
                atom2 = &res2->atoms[0];
            }
        }
        if (!atom2) atom2 = atom1;

        coot::util::water_coordination_t wc;
        coot::util::contact_atoms_info_t cai(atom1, coot::util::contact_atom_t(atom2, atom1));
        wc.atom_contacts.push_back(std::move(cai));

        float orig_x = atom2->pos.x;
        float orig_y = atom2->pos.y;
        float orig_z = atom2->pos.z;

        atom2->pos.x = orig_x;
        atom2->pos.y = orig_y;
        atom2->pos.z = orig_z;

        wc.atom_contacts[0][0].mat[3][0] = 10.0;
        wc.atom_contacts[0][0].mat[3][1] = 0.0;
        wc.atom_contacts[0][0].mat[3][2] = 0.0;
        wc.atom_contacts[0][0].mat[3][3] = 1.0;

        wc.transform_atom_gemmi(0, 0);

        EXPECT_NEAR(atom2->pos.x, 13.401, 1e-4);
        EXPECT_NEAR(atom2->pos.y, 3.78, 1e-4);
        EXPECT_NEAR(atom2->pos.z, 74.703, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}