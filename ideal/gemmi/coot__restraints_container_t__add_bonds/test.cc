#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <vector>
#include <set>
#include <map>
#include <string>
#include <algorithm>

TEST(OracleTest, add_bonds) {
    // Case 1: ALA 21 (original oracle case)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        coot::protein_geometry geom;
        geom.init_standard();

        // Find residue A/21
        gemmi::Model& model = st.models[0];
        gemmi::Residue* residue_p = nullptr;

        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 21) {
                        residue_p = &res;
                        break;
                    }
                }
            }
            if (residue_p) break;
        }

        EXPECT_NE(residue_p, nullptr);
        EXPECT_EQ(residue_p->seqid.num.value, 21);

        // Collect atoms from this residue
        std::vector<gemmi::Atom*> res_atoms;
        for (gemmi::Atom& atom : residue_p->atoms) {
            res_atoms.push_back(&atom);
        }
        EXPECT_EQ(res_atoms.size(), 10u);

        std::string pdb_resname(residue_p->name);
        int idr = geom.get_monomer_restraints_index(pdb_resname, -1, false);
        EXPECT_EQ(idr, 0);

        std::vector<coot::simple_restraint> restraints_vec;
        std::vector<std::set<int>> bonded_atom_indices;
        std::map<gemmi::Atom*, coot::hb_t> H_atom_parent_energy_type_atom_map;

        int n_restraints_before = restraints_vec.size();
        EXPECT_EQ(n_restraints_before, 0);

        int n_bond_restr = coot::add_bonds_gemmi(restraints_vec, bonded_atom_indices, H_atom_parent_energy_type_atom_map,
                                                  idr, res_atoms, residue_p, geom);

        int n_restraints_after = restraints_vec.size();

        EXPECT_EQ(n_bond_restr, 9);
        EXPECT_EQ(n_restraints_after, 9);
        EXPECT_EQ(n_restraints_after - n_restraints_before, 9);
    }

    // Case 2: GLY 8 - different residue (no sidechain, fewer atoms/bonds)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        coot::protein_geometry geom;
        geom.init_standard();

        // Find residue A/8
        gemmi::Model& model = st.models[0];
        gemmi::Residue* residue_p = nullptr;

        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 8) {
                        residue_p = &res;
                        break;
                    }
                }
            }
            if (residue_p) break;
        }

        EXPECT_NE(residue_p, nullptr);
        EXPECT_EQ(residue_p->seqid.num.value, 8);

        // Collect atoms from this residue
        std::vector<gemmi::Atom*> res_atoms;
        for (gemmi::Atom& atom : residue_p->atoms) {
            res_atoms.push_back(&atom);
        }
        EXPECT_GT(res_atoms.size(), 0u);

        std::string pdb_resname(residue_p->name);
        int idr = geom.get_monomer_restraints_index(pdb_resname, -1, false);
        EXPECT_GE(idr, 0);

        std::vector<coot::simple_restraint> restraints_vec;
        std::vector<std::set<int>> bonded_atom_indices;
        std::map<gemmi::Atom*, coot::hb_t> H_atom_parent_energy_type_atom_map;

        int n_restraints_before = restraints_vec.size();
        EXPECT_EQ(n_restraints_before, 0);

        int n_bond_restr = coot::add_bonds_gemmi(restraints_vec, bonded_atom_indices, H_atom_parent_energy_type_atom_map,
                                                  idr, res_atoms, residue_p, geom);

        int n_restraints_after = restraints_vec.size();

        EXPECT_GT(n_bond_restr, 0);
        EXPECT_GT(n_restraints_after, 0);
        EXPECT_EQ(n_restraints_after - n_restraints_before, n_bond_restr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}