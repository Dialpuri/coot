#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static std::string trim(const std::string& s) {
    int i = 0; int j = static_cast<int>(s.size()) - 1;
    while (i <= j && s[i] == ' ') i++;
    while (j >= i && s[j] == ' ') j--;
    return s.substr(i, j - i + 1);
}

static gemmi::Atom* find_atom(gemmi::Residue& res, const std::string& name) {
    std::string trimmed_name = trim(name);
    for (gemmi::Atom& a : res.atoms) {
        if (trim(a.name) == trimmed_name)
            return &a;
    }
    return nullptr;
}

TEST(OracleTest, add_main_chain_hydrogens) {
    {
        // Case 1: ALA (non-Gly, non-Pro) — should call add_main_chain_HA and add_main_chain_H
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];

        gemmi::Residue* ala_res = nullptr;
        gemmi::Residue* prev_res = nullptr;
        for (size_t i = 0; i < chain.residues.size(); i++) {
            gemmi::Residue& r = chain.residues[i];
            if (r.name == "ALA") {
                ala_res = &r;
                if (i > 0) {
                    prev_res = &chain.residues[i - 1];
                }
                break;
            }
        }

        EXPECT_NE(ala_res, nullptr);
        EXPECT_EQ(ala_res->name, "ALA");
        EXPECT_EQ(ala_res->seqid.num.value, 21);
        EXPECT_NE(prev_res, nullptr);
        EXPECT_EQ(prev_res->name, "LYS");

        gemmi::Atom* ha_atom = find_atom(*ala_res, " HA ");
        double ha_x_before = 0, ha_y_before = 0, ha_z_before = 0;
        if (ha_atom) {
            ha_x_before = ha_atom->pos.x;
            ha_y_before = ha_atom->pos.y;
            ha_z_before = ha_atom->pos.z;
        }
        EXPECT_NEAR(ha_x_before, 23.908, 1e-4);
        EXPECT_NEAR(ha_y_before, 6.621, 1e-4);
        EXPECT_NEAR(ha_z_before, 75.009, 1e-4);

        gemmi::Atom* hn_atom = find_atom(*ala_res, " H  ");
        double hn_x_before = 0, hn_y_before = 0, hn_z_before = 0;
        if (hn_atom) {
            hn_x_before = hn_atom->pos.x;
            hn_y_before = hn_atom->pos.y;
            hn_z_before = hn_atom->pos.z;
        }
        EXPECT_NEAR(hn_x_before, 25.12, 1e-4);
        EXPECT_NEAR(hn_y_before, 5.495, 1e-4);
        EXPECT_NEAR(hn_z_before, 72.747, 1e-4);

        coot::reduce::add_main_chain_hydrogens_gemmi(*ala_res, prev_res, false);

        gemmi::Atom* ha_atom_after = find_atom(*ala_res, " HA ");
        double ha_x_after = 0, ha_y_after = 0, ha_z_after = 0;
        if (ha_atom_after) {
            ha_x_after = ha_atom_after->pos.x;
            ha_y_after = ha_atom_after->pos.y;
            ha_z_after = ha_atom_after->pos.z;
        }
        EXPECT_NEAR(ha_x_after, 23.8997, 1e-4);
        EXPECT_NEAR(ha_y_after, 6.60619, 1e-4);
        EXPECT_NEAR(ha_z_after, 74.9858, 1e-4);

        gemmi::Atom* hn_atom_after = find_atom(*ala_res, " H  ");
        double hn_x_after = 0, hn_y_after = 0, hn_z_after = 0;
        if (hn_atom_after) {
            hn_x_after = hn_atom_after->pos.x;
            hn_y_after = hn_atom_after->pos.y;
            hn_z_after = hn_atom_after->pos.z;
        }
        EXPECT_NEAR(hn_x_after, 25.061, 1e-4);
        EXPECT_NEAR(hn_y_after, 5.48895, 1e-4);
        EXPECT_NEAR(hn_z_after, 72.8198, 1e-4);
    }

    {
        // Case 2: GLY (is_gly=true) — should call add_2_sp3_hydrogens for HA2/HA3 + add_main_chain_H
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];

        gemmi::Residue* gly_res = nullptr;
        gemmi::Residue* gly_prev = nullptr;
        for (size_t i = 0; i < chain.residues.size(); i++) {
            gemmi::Residue& r = chain.residues[i];
            if (r.name == "GLY") {
                gly_res = &r;
                if (i > 0) {
                    gly_prev = &chain.residues[i - 1];
                }
                break;
            }
        }

        EXPECT_NE(gly_res, nullptr);
        EXPECT_EQ(gly_res->name, "GLY");
        EXPECT_EQ(gly_res->seqid.num.value, 11);
        EXPECT_NE(gly_prev, nullptr);
        EXPECT_EQ(gly_prev->name, "ILE");

        gemmi::Atom* ha2_atom = find_atom(*gly_res, " HA2");
        double ha2_x_before = 0, ha2_y_before = 0, ha2_z_before = 0;
        if (ha2_atom) {
            ha2_x_before = ha2_atom->pos.x;
            ha2_y_before = ha2_atom->pos.y;
            ha2_z_before = ha2_atom->pos.z;
        }
        EXPECT_NEAR(ha2_x_before, 25.876, 1e-4);
        EXPECT_NEAR(ha2_y_before, 14.76, 1e-4);
        EXPECT_NEAR(ha2_z_before, 64.703, 1e-4);

        gemmi::Atom* ha3_atom = find_atom(*gly_res, " HA3");
        double ha3_x_before = 0, ha3_y_before = 0, ha3_z_before = 0;
        if (ha3_atom) {
            ha3_x_before = ha3_atom->pos.x;
            ha3_y_before = ha3_atom->pos.y;
            ha3_z_before = ha3_atom->pos.z;
        }
        EXPECT_NEAR(ha3_x_before, 25.943, 1e-4);
        EXPECT_NEAR(ha3_y_before, 13.322, 1e-4);
        EXPECT_NEAR(ha3_z_before, 64.013, 1e-4);

        coot::reduce::add_main_chain_hydrogens_gemmi(*gly_res, gly_prev, true);

        gemmi::Atom* ha2_atom_after = find_atom(*gly_res, " HA2");
        double ha2_x_after = 0, ha2_y_after = 0, ha2_z_after = 0;
        if (ha2_atom_after) {
            ha2_x_after = ha2_atom_after->pos.x;
            ha2_y_after = ha2_atom_after->pos.y;
            ha2_z_after = ha2_atom_after->pos.z;
        }
        EXPECT_NEAR(ha2_x_after, 25.8719, 1e-4);
        EXPECT_NEAR(ha2_y_after, 14.7442, 1e-4);
        EXPECT_NEAR(ha2_z_after, 64.693, 1e-4);

        gemmi::Atom* ha3_atom_after = find_atom(*gly_res, " HA3");
        double ha3_x_after = 0, ha3_y_after = 0, ha3_z_after = 0;
        if (ha3_atom_after) {
            ha3_x_after = ha3_atom_after->pos.x;
            ha3_y_after = ha3_atom_after->pos.y;
            ha3_z_after = ha3_atom_after->pos.z;
        }
        EXPECT_NEAR(ha3_x_after, 25.9376, 1e-4);
        EXPECT_NEAR(ha3_y_after, 13.3395, 1e-4);
        EXPECT_NEAR(ha3_z_after, 64.0192, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
