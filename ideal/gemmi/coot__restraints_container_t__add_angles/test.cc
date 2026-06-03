#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <vector>
#include <set>
#include <string>

TEST(OracleTest, add_angles) {
    // Case 1: LYS at chain A, residue 6
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        coot::protein_geometry geom;
        geom.init_standard();

        // Get a residue — LYS at chain A, residue 6
        gemmi::Model& model = st.models[0];
        gemmi::Residue* residue_p = nullptr;

        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 6) {
                        residue_p = &res;
                        break;
                    }
                }
            }
            if (residue_p) break;
        }

        ASSERT_NE(residue_p, nullptr);

        // Get atoms from the residue
        std::vector<gemmi::Atom*> res_selection;
        for (gemmi::Atom& atom : residue_p->atoms) {
            res_selection.push_back(&atom);
        }

        // Get idr (index dictionary residue)
        std::string pdb_resname(residue_p->name);
        int idr = geom.get_monomer_restraints_index(pdb_resname, -1, false);

        std::vector<coot::simple_restraint> restraints_vec;
        std::vector<std::set<int>> bonded_atom_indices;
        std::set<int> fixed_atom_indices;

        int n_before = restraints_vec.size();

        // Call add_angles directly
        int n_added = coot::add_angles_gemmi(restraints_vec, bonded_atom_indices, fixed_atom_indices, idr, res_selection, residue_p, geom);

        int n_after = restraints_vec.size();

        EXPECT_EQ(n_added, 38);
        EXPECT_EQ(n_before, 0);
        EXPECT_EQ(n_after, 38);
    }

    // Case 2: MET at chain A, residue 1
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        coot::protein_geometry geom;
        geom.init_standard();

        gemmi::Model& model = st.models[0];

        // First add LYS restraints (to match oracle state where tyr_before=38)
        gemmi::Residue* residue_lys = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 6) {
                        residue_lys = &res;
                        break;
                    }
                }
            }
            if (residue_lys) break;
        }
        ASSERT_NE(residue_lys, nullptr);

        std::vector<gemmi::Atom*> lys_selection;
        for (gemmi::Atom& atom : residue_lys->atoms) {
            lys_selection.push_back(&atom);
        }

        std::string lys_name(residue_lys->name);
        int idr_lys = geom.get_monomer_restraints_index(lys_name, -1, false);

        std::vector<coot::simple_restraint> restraints_vec;
        std::vector<std::set<int>> bonded_atom_indices;
        std::set<int> fixed_atom_indices;

        // Add LYS angles first (to get restraints_before = 38)
        coot::add_angles_gemmi(restraints_vec, bonded_atom_indices, fixed_atom_indices, idr_lys, lys_selection, residue_lys, geom);

        // Now add MET angles
        gemmi::Residue* residue_met = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 1) {
                        residue_met = &res;
                        break;
                    }
                }
            }
            if (residue_met) break;
        }
        ASSERT_NE(residue_met, nullptr);

        std::string met_name(residue_met->name);
        int idr_met = geom.get_monomer_restraints_index(met_name, -1, false);

        std::vector<gemmi::Atom*> met_selection;
        for (gemmi::Atom& atom : residue_met->atoms) {
            met_selection.push_back(&atom);
        }

        int n_met_before = restraints_vec.size();
        int n_met_added = coot::add_angles_gemmi(restraints_vec, bonded_atom_indices, fixed_atom_indices, idr_met, met_selection, residue_met, geom);
        int n_met_after = restraints_vec.size();

        EXPECT_EQ(n_met_added, 32);
        EXPECT_EQ(n_met_before, 38);
        EXPECT_EQ(n_met_after, 70);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}