#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a residue by chain name and sequence number in a gemmi model
static gemmi::Residue* find_residue(gemmi::Model& model, std::string chain_name, int seq_num) {
   for (gemmi::Chain& ch : model.chains) {
      if (ch.name == chain_name) {
         for (gemmi::Residue& res : ch.residues) {
            if (res.seqid.num.value == seq_num) {
               return &res;
            }
         }
      }
   }
   return nullptr;
}

TEST(OracleTest, torsionable_bonds_monomer_internal) {
    // Load structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];

    // Initialize protein geometry
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: Valid residue with side-chain torsions — MET 1
    {
        gemmi::Residue *res = find_residue(model, "A", 1);
        ASSERT_NE(res, nullptr);

        std::vector<gemmi::Atom*> atom_sel;
        for (gemmi::Atom& atom : res->atoms) {
            atom_sel.push_back(&atom);
        }
        int n_selected = static_cast<int>(atom_sel.size());

        ASSERT_EQ(res->name, "MET");
        ASSERT_EQ(res->seqid.num.value, 1);
        EXPECT_EQ(n_selected, 19);

        std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> result =
            coot::torsionable_bonds_monomer_internal_gemmi(res, atom_sel, n_selected, false, &geom);

        EXPECT_EQ(result.size(), 6u);
        if (result.size() >= 1) EXPECT_EQ(coot::pad_atom_name(result[0].first->name, result[0].first->element.name()), " CA ");
        if (result.size() >= 1) EXPECT_EQ(coot::pad_atom_name(result[0].second->name, result[0].second->element.name()), " CB ");
        if (result.size() >= 2) EXPECT_EQ(coot::pad_atom_name(result[1].first->name, result[1].first->element.name()), " CB ");
        if (result.size() >= 2) EXPECT_EQ(coot::pad_atom_name(result[1].second->name, result[1].second->element.name()), " CG ");
        if (result.size() >= 3) EXPECT_EQ(coot::pad_atom_name(result[2].first->name, result[2].first->element.name()), " CG ");
        if (result.size() >= 3) EXPECT_EQ(coot::pad_atom_name(result[2].second->name, result[2].second->element.name()), " SD ");
        if (result.size() >= 4) EXPECT_EQ(coot::pad_atom_name(result[3].first->name, result[3].first->element.name()), " SD ");
        if (result.size() >= 4) EXPECT_EQ(coot::pad_atom_name(result[3].second->name, result[3].second->element.name()), " CE ");
        if (result.size() >= 5) EXPECT_EQ(coot::pad_atom_name(result[4].first->name, result[4].first->element.name()), " CA ");
        if (result.size() >= 5) EXPECT_EQ(coot::pad_atom_name(result[4].second->name, result[4].second->element.name()), " N  ");
        if (result.size() >= 6) EXPECT_EQ(coot::pad_atom_name(result[5].first->name, result[5].first->element.name()), " C  ");
        if (result.size() >= 6) EXPECT_EQ(coot::pad_atom_name(result[5].second->name, result[5].second->element.name()), " CA ");
    }

    // Case 2: Residue with no side-chain torsions — GLY 11 (only backbone)
    {
        gemmi::Residue *res = find_residue(model, "A", 11);
        ASSERT_NE(res, nullptr);

        std::vector<gemmi::Atom*> atom_sel;
        for (gemmi::Atom& atom : res->atoms) {
            atom_sel.push_back(&atom);
        }
        int n_selected = static_cast<int>(atom_sel.size());

        ASSERT_EQ(res->name, "GLY");
        ASSERT_EQ(res->seqid.num.value, 11);
        EXPECT_EQ(n_selected, 7);

        std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> result =
            coot::torsionable_bonds_monomer_internal_gemmi(res, atom_sel, n_selected, false, &geom);

        EXPECT_EQ(result.size(), 2u);
        if (result.size() >= 1) EXPECT_EQ(coot::pad_atom_name(result[0].first->name, result[0].first->element.name()), " CA ");
        if (result.size() >= 1) EXPECT_EQ(coot::pad_atom_name(result[0].second->name, result[0].second->element.name()), " N  ");
        if (result.size() >= 2) EXPECT_EQ(coot::pad_atom_name(result[1].first->name, result[1].first->element.name()), " C  ");
        if (result.size() >= 2) EXPECT_EQ(coot::pad_atom_name(result[1].second->name, result[1].second->element.name()), " CA ");
    }

    // Case 3: Empty atom selection — should return empty
    {
        gemmi::Residue *res = find_residue(model, "A", 1);
        ASSERT_NE(res, nullptr);

        std::vector<gemmi::Atom*> empty_sel;
        int n_selected = 0;

        ASSERT_EQ(res->name, "MET");
        ASSERT_EQ(res->seqid.num.value, 1);

        std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> result =
            coot::torsionable_bonds_monomer_internal_gemmi(res, empty_sel, n_selected, false, &geom);

        EXPECT_EQ(result.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}