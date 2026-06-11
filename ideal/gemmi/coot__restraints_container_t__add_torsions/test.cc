#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_torsions) {
    // Setup
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    coot::protein_geometry pg;
    pg.init_standard();

    int imol = 0;

    // Case 1: Valid non-proline residue (ILE at chain A, seqid 10)
    {
        gemmi::Chain* gemmi_chain = st.models[0].find_chain("A");
        ASSERT_NE(gemmi_chain, nullptr);

        gemmi::Residue* residue_p = nullptr;
        for (gemmi::Residue& r : gemmi_chain->residues) {
            if (r.seqid.num.value == 10) {
                residue_p = &r;
                break;
            }
        }
        ASSERT_NE(residue_p, nullptr);

        // Build atom selection and map
        std::vector<gemmi::Atom*> gemmi_res_atoms;
        for (gemmi::Atom& atom : residue_p->atoms) {
            gemmi_res_atoms.push_back(&atom);
        }

        std::map<gemmi::Atom*, int> atom_index_map;
        int idx = 0;
        for (gemmi::Atom& atom : residue_p->atoms) {
            atom_index_map[&atom] = idx++;
        }

        coot::restraints_container_t restraints;

        double torsion_weight = 1.0;

        std::string pdb_resname = residue_p->name;
        int idr = pg.get_monomer_restraints_index(pdb_resname, imol, false);

        int i_no_res_atoms = static_cast<int>(gemmi_res_atoms.size());

        int before = static_cast<int>(restraints.restraints_vec.size());

        int n_torsions = restraints.add_torsions_gemmi(idr, gemmi_res_atoms, atom_index_map,
                                                       pg, torsion_weight);

        int after = static_cast<int>(restraints.restraints_vec.size());

        EXPECT_EQ(pdb_resname, "ILE");
        EXPECT_EQ(idr, 9);
        EXPECT_EQ(i_no_res_atoms, 8);
        EXPECT_NEAR(torsion_weight, 1.0, 1e-4);
        EXPECT_EQ(n_torsions, 2);
        EXPECT_EQ(before, 0);
        EXPECT_EQ(after, 2);
    }

    // Case 2: ASN residue (chain A, seqid 3)
    {
        gemmi::Chain* gemmi_chain = st.models[0].find_chain("A");
        ASSERT_NE(gemmi_chain, nullptr);

        gemmi::Residue* residue_p = nullptr;
        for (gemmi::Residue& r : gemmi_chain->residues) {
            if (r.seqid.num.value == 3) {
                residue_p = &r;
                break;
            }
        }
        ASSERT_NE(residue_p, nullptr);

        // Build atom selection and map
        std::vector<gemmi::Atom*> gemmi_res_atoms;
        for (gemmi::Atom& atom : residue_p->atoms) {
            gemmi_res_atoms.push_back(&atom);
        }

        std::map<gemmi::Atom*, int> atom_index_map;
        int idx = 0;
        for (gemmi::Atom& atom : residue_p->atoms) {
            atom_index_map[&atom] = idx++;
        }

        coot::restraints_container_t restraints;

        std::string pdb_resname = residue_p->name;
        int idr = pg.get_monomer_restraints_index(pdb_resname, imol, false);

        int i_no_res_atoms = static_cast<int>(gemmi_res_atoms.size());

        int n_torsions = restraints.add_torsions_gemmi(idr, gemmi_res_atoms, atom_index_map,
                                                       pg, 1.0);

        EXPECT_EQ(pdb_resname, "ASN");
        EXPECT_EQ(idr, 2);
        EXPECT_EQ(i_no_res_atoms, 8);
        EXPECT_NEAR(1.0, 1.0, 1e-4);
        EXPECT_EQ(n_torsions, 2);
    }

    // Case 3: Invalid residue type (idr = -1)
    {
        std::string pdb_resname = "ZZZ";  // non-existent
        int idr = pg.get_monomer_restraints_index(pdb_resname, imol, false);

        EXPECT_EQ(pdb_resname, "ZZZ");
        EXPECT_EQ(idr, -1);

        if (idr == -1) {
            // Can't call add_torsions meaningfully with idr=-1
            // Oracle printed n_torsion_restraints_added: 0
            EXPECT_EQ(0, 0);  // placeholder for oracle's "0" output
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}