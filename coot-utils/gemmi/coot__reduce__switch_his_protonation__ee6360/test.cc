#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper to find a residue by chain name and sequence number
gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seq_num) {
    for (auto& chain : model.chains) {
        if (chain.name == chain_id) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

// Helper to find an atom by name in a residue (non-const)
gemmi::Atom* find_atom(gemmi::Residue& res, const std::string& name) {
    for (auto& a : res.atoms) {
        if (a.name == name) {
            return &a;
        }
    }
    return nullptr;
}

TEST(OracleTest, switch_his_protonation) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Case 1: Valid HIS residue with HD1 -> HE2 switch
    {
        gemmi::Residue* res = find_residue(st.models[0], "A", 60);
        ASSERT_NE(res, nullptr);

        // Check required atoms exist
        gemmi::Atom* cg  = find_atom(*res, "CG");
        gemmi::Atom* nd1 = find_atom(*res, "ND1");
        gemmi::Atom* ce1 = find_atom(*res, "CE1");

        EXPECT_TRUE(cg);
        EXPECT_TRUE(nd1);
        EXPECT_TRUE(ce1);

        // Create an HD1 hydrogen atom to serve as current_H_atom
        gemmi::Atom h_atom;
        h_atom.name = "HD1";
        h_atom.element = gemmi::El::H;
        h_atom.pos = gemmi::Position(0.0, 0.0, 0.0);
        h_atom.occ = 1.0;
        h_atom.b_iso = 10.0;
        h_atom.altloc = '\0';
        res->atoms.push_back(std::move(h_atom));

        gemmi::Atom* h_atom_ptr = &res->atoms.back();
        EXPECT_STREQ(h_atom_ptr->name.c_str(), "HD1");

        // Call the function
        EXPECT_NO_THROW(coot::reduce::switch_his_protonation_gemmi(res, h_atom_ptr, 0.93));

        EXPECT_STREQ(h_atom_ptr->name.c_str(), "HE2");
        EXPECT_NEAR(h_atom_ptr->pos.x, 31.5227, 1e-4);
        EXPECT_NEAR(h_atom_ptr->pos.y, -8.11556, 1e-4);
        EXPECT_NEAR(h_atom_ptr->pos.z, 53.9955, 1e-4);
    }

    // Case 2: Null current_H_atom (guarded path)
    {
        gemmi::Residue* res = find_residue(st.models[0], "A", 60);

        EXPECT_NO_THROW(coot::reduce::switch_his_protonation_gemmi(res, nullptr, 0.93));
    }

    // Case 3: HE2 -> HD1 switch
    {
        gemmi::Residue* res = find_residue(st.models[0], "A", 71);
        ASSERT_NE(res, nullptr);

        gemmi::Atom* ce1 = find_atom(*res, "CE1");
        gemmi::Atom* ne2 = find_atom(*res, "NE2");
        gemmi::Atom* cd2 = find_atom(*res, "CD2");

        EXPECT_TRUE(ce1);
        EXPECT_TRUE(ne2);
        EXPECT_TRUE(cd2);

        gemmi::Atom h_atom;
        h_atom.name = "HE2";
        h_atom.element = gemmi::El::H;
        h_atom.pos = gemmi::Position(0.0, 0.0, 0.0);
        h_atom.occ = 1.0;
        h_atom.b_iso = 10.0;
        h_atom.altloc = '\0';
        res->atoms.push_back(std::move(h_atom));

        gemmi::Atom* h_atom_ptr = &res->atoms.back();
        EXPECT_STREQ(h_atom_ptr->name.c_str(), "HE2");

        EXPECT_NO_THROW(coot::reduce::switch_his_protonation_gemmi(res, h_atom_ptr, 0.93));

        EXPECT_STREQ(h_atom_ptr->name.c_str(), "HD1");
        EXPECT_NEAR(h_atom_ptr->pos.x, 5.65211, 1e-4);
        EXPECT_NEAR(h_atom_ptr->pos.y, 0.639522, 1e-4);
        EXPECT_NEAR(h_atom_ptr->pos.z, 64.6046, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}