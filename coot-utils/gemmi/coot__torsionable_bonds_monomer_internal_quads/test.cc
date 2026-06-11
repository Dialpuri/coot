#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "coot-utils/coot-coord-utils.hh"

#include <vector>
#include <string>

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

static std::string pad_name(const gemmi::Atom* atom) {
    return coot::pad_atom_name(atom->name, atom->element.name());
}

TEST(OracleTest, torsionable_bonds_monomer_internal_quads) {
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: MET residue (first residue in chain A)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        gemmi::Model &model = st.models[0];
        gemmi::Chain &chain = model.chains[0];

        gemmi::Residue &res = chain.residues[0];
        EXPECT_NE(&res, nullptr);

        std::vector<gemmi::Atom*> atom_table;
        for (gemmi::Atom &a : res.atoms) {
            atom_table.push_back(&a);
        }
        int n_atoms = static_cast<int>(atom_table.size());

        EXPECT_EQ(std::string(res.name), "MET");
        EXPECT_EQ(std::string(chain.name), "A");
        EXPECT_EQ(res.seqid.num.value, 1);
        EXPECT_EQ(n_atoms, 19);

        std::vector<coot::torsion_atom_quad_gemmi> quads =
            coot::torsionable_bonds_monomer_internal_quads_gemmi(res, atom_table, false, &geom);

        EXPECT_EQ(quads.size(), 6u);

        // quad[0]: chi1
        EXPECT_EQ(quads[0].name, "chi1");
        EXPECT_EQ(quads[0].residue_name, "MET");
        EXPECT_NEAR(quads[0].angle, -60.0, 1e-4);
        EXPECT_NEAR(quads[0].angle_esd, 10.0, 1e-4);
        EXPECT_EQ(quads[0].period, 3);
        EXPECT_NE(quads[0].atom_1, nullptr);
        EXPECT_EQ(pad_name(quads[0].atom_1), " N  ");
        EXPECT_EQ(pad_name(quads[0].atom_2), " CA ");
        EXPECT_EQ(pad_name(quads[0].atom_3), " CB ");
        EXPECT_EQ(pad_name(quads[0].atom_4), " CG ");

        // quad[1]: chi2
        EXPECT_EQ(quads[1].name, "chi2");
        EXPECT_EQ(quads[1].residue_name, "MET");
        EXPECT_NEAR(quads[1].angle, 180.0, 1e-4);
        EXPECT_NEAR(quads[1].angle_esd, 10.0, 1e-4);
        EXPECT_EQ(quads[1].period, 3);
        EXPECT_EQ(pad_name(quads[1].atom_1), " CA ");
        EXPECT_EQ(pad_name(quads[1].atom_2), " CB ");
        EXPECT_EQ(pad_name(quads[1].atom_3), " CG ");
        EXPECT_EQ(pad_name(quads[1].atom_4), " SD ");

        // quad[2]: chi3
        EXPECT_EQ(quads[2].name, "chi3");
        EXPECT_EQ(quads[2].residue_name, "MET");
        EXPECT_NEAR(quads[2].angle, 180.0, 1e-4);
        EXPECT_NEAR(quads[2].angle_esd, 10.0, 1e-4);
        EXPECT_EQ(quads[2].period, 3);
        EXPECT_EQ(pad_name(quads[2].atom_1), " CB ");
        EXPECT_EQ(pad_name(quads[2].atom_2), " CG ");
        EXPECT_EQ(pad_name(quads[2].atom_3), " SD ");
        EXPECT_EQ(pad_name(quads[2].atom_4), " CE ");

        // quad[3]: hh1
        EXPECT_EQ(quads[3].name, "hh1");
        EXPECT_EQ(quads[3].residue_name, "MET");
        EXPECT_NEAR(quads[3].angle, 60.0, 1e-4);
        EXPECT_NEAR(quads[3].angle_esd, 10.0, 1e-4);
        EXPECT_EQ(quads[3].period, 3);
        EXPECT_EQ(pad_name(quads[3].atom_1), " CG ");
        EXPECT_EQ(pad_name(quads[3].atom_2), " SD ");
        EXPECT_EQ(pad_name(quads[3].atom_3), " CE ");
        EXPECT_EQ(pad_name(quads[3].atom_4), " HE3");

        // quad[4]: sp3_sp3_1
        EXPECT_EQ(quads[4].name, "sp3_sp3_1");
        EXPECT_EQ(quads[4].residue_name, "MET");
        EXPECT_NEAR(quads[4].angle, 180.0, 1e-4);
        EXPECT_NEAR(quads[4].angle_esd, 10.0, 1e-4);
        EXPECT_EQ(quads[4].period, 3);
        EXPECT_EQ(pad_name(quads[4].atom_1), " C  ");
        EXPECT_EQ(pad_name(quads[4].atom_2), " CA ");
        EXPECT_EQ(pad_name(quads[4].atom_3), " N  ");
        EXPECT_EQ(pad_name(quads[4].atom_4), " H  ");

        // quad[5]: sp2_sp3_1
        EXPECT_EQ(quads[5].name, "sp2_sp3_1");
        EXPECT_EQ(quads[5].residue_name, "MET");
        EXPECT_NEAR(quads[5].angle, 0.0, 1e-4);
        EXPECT_NEAR(quads[5].angle_esd, 10.0, 1e-4);
        EXPECT_EQ(quads[5].period, 6);
        EXPECT_EQ(pad_name(quads[5].atom_1), " O  ");
        EXPECT_EQ(pad_name(quads[5].atom_2), " C  ");
        EXPECT_EQ(pad_name(quads[5].atom_3), " CA ");
        EXPECT_EQ(pad_name(quads[5].atom_4), " N  ");
    }

    // Case 2: HOH residue (chain B) - should have no torsions
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        gemmi::Model &model = st.models[0];
        gemmi::Chain &chain = model.chains[1];
        EXPECT_NE(&chain, nullptr);

        gemmi::Residue &res2 = chain.residues[0];
        EXPECT_NE(&res2, nullptr);

        std::vector<gemmi::Atom*> atom_table2;
        for (gemmi::Atom &a : res2.atoms) {
            atom_table2.push_back(&a);
        }
        int n_atoms2 = static_cast<int>(atom_table2.size());

        EXPECT_EQ(std::string(res2.name), "HOH");
        EXPECT_EQ(std::string(chain.name), "B");
        EXPECT_EQ(res2.seqid.num.value, 1);
        EXPECT_EQ(n_atoms2, 1);

        std::vector<coot::torsion_atom_quad_gemmi> quads2 =
            coot::torsionable_bonds_monomer_internal_quads_gemmi(res2, atom_table2, false, &geom);

        EXPECT_EQ(quads2.size(), 0u);
    }

    // Case 3: HIS residue with include_pyranose_ring_torsions=true
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        gemmi::Model &model = st.models[0];
        gemmi::Chain &chain = model.chains[0];

        gemmi::Residue *res3 = nullptr;
        for (gemmi::Residue &r : chain.residues) {
            if (r.name == "HIS") {
                res3 = &r;
                break;
            }
        }
        EXPECT_NE(res3, nullptr);

        std::vector<gemmi::Atom*> atom_table3;
        for (gemmi::Atom &a : res3->atoms) {
            atom_table3.push_back(&a);
        }
        int n_atoms3 = static_cast<int>(atom_table3.size());

        EXPECT_EQ(std::string(res3->name), "HIS");
        EXPECT_EQ(std::string(chain.name), "A");
        EXPECT_EQ(res3->seqid.num.value, 60);
        EXPECT_EQ(n_atoms3, 18);

        std::vector<coot::torsion_atom_quad_gemmi> quads3 =
            coot::torsionable_bonds_monomer_internal_quads_gemmi(*res3, atom_table3, true, &geom);

        EXPECT_EQ(quads3.size(), 4u);

        // quad[0]: chi1
        EXPECT_EQ(quads3[0].name, "chi1");
        EXPECT_EQ(quads3[0].residue_name, "HIS");
        EXPECT_NEAR(quads3[0].angle, -60.0, 1e-4);
        EXPECT_EQ(pad_name(quads3[0].atom_1), " N  ");
        EXPECT_EQ(pad_name(quads3[0].atom_2), " CA ");
        EXPECT_EQ(pad_name(quads3[0].atom_3), " CB ");
        EXPECT_EQ(pad_name(quads3[0].atom_4), " CG ");

        // quad[1]: chi2
        EXPECT_EQ(quads3[1].name, "chi2");
        EXPECT_EQ(quads3[1].residue_name, "HIS");
        EXPECT_NEAR(quads3[1].angle, 90.0, 1e-4);
        EXPECT_EQ(pad_name(quads3[1].atom_1), " CA ");
        EXPECT_EQ(pad_name(quads3[1].atom_2), " CB ");
        EXPECT_EQ(pad_name(quads3[1].atom_3), " CG ");
        EXPECT_EQ(pad_name(quads3[1].atom_4), " CD2");

        // quad[2]: sp3_sp3_1
        EXPECT_EQ(quads3[2].name, "sp3_sp3_1");
        EXPECT_EQ(quads3[2].residue_name, "HIS");
        EXPECT_NEAR(quads3[2].angle, 180.0, 1e-4);
        EXPECT_EQ(pad_name(quads3[2].atom_1), " C  ");
        EXPECT_EQ(pad_name(quads3[2].atom_2), " CA ");
        EXPECT_EQ(pad_name(quads3[2].atom_3), " N  ");
        EXPECT_EQ(pad_name(quads3[2].atom_4), " H  ");

        // quad[3]: sp2_sp3_1
        EXPECT_EQ(quads3[3].name, "sp2_sp3_1");
        EXPECT_EQ(quads3[3].residue_name, "HIS");
        EXPECT_NEAR(quads3[3].angle, 0.0, 1e-4);
        EXPECT_EQ(pad_name(quads3[3].atom_1), " O  ");
        EXPECT_EQ(pad_name(quads3[3].atom_2), " C  ");
        EXPECT_EQ(pad_name(quads3[3].atom_3), " CA ");
        EXPECT_EQ(pad_name(quads3[3].atom_4), " N  ");
    }
}
