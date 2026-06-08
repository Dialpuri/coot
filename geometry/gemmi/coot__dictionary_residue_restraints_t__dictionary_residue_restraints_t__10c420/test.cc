#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, dictionary_residue_restraints_t) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model &model = st.models[0];

    // Case 1: ASN at chain A (GetResidue(2) → seqnum 3)
    {
        gemmi::Chain *chain = model.find_chain("A");
        ASSERT_NE(chain, nullptr);
        gemmi::Residue *res = &chain->residues[2];
        ASSERT_NE(res, nullptr);
        EXPECT_EQ(res->seqid.num.value, 3);
        EXPECT_EQ(res->name, "ASN");
        EXPECT_EQ(res->atoms.size(), 8u);

        coot::dictionary_residue_restraints_t dict_res =
            coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t_gemmi(res);

        EXPECT_EQ(dict_res.residue_info.comp_id, "ASN");
        EXPECT_EQ(dict_res.atom_info.size(), 8u);
        EXPECT_EQ(dict_res.bond_restraint.size(), 7u);
        EXPECT_EQ(dict_res.angle_restraint.size(), 8u);
        EXPECT_EQ(dict_res.chiral_restraint.size(), 0u);
        EXPECT_EQ(dict_res.plane_restraint.size(), 0u);
    }

    // Case 2: PHE at chain A (GetResidue(3) → seqnum 4)
    {
        gemmi::Chain *chain = model.find_chain("A");
        gemmi::Residue *res = &chain->residues[3];
        ASSERT_NE(res, nullptr);
        EXPECT_EQ(res->seqid.num.value, 4);
        EXPECT_EQ(res->name, "PHE");
        EXPECT_EQ(res->atoms.size(), 11u);

        coot::dictionary_residue_restraints_t dict_res =
            coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t_gemmi(res);

        EXPECT_EQ(dict_res.residue_info.comp_id, "PHE");
        EXPECT_EQ(dict_res.atom_info.size(), 11u);
        EXPECT_EQ(dict_res.bond_restraint.size(), 11u);
        EXPECT_EQ(dict_res.angle_restraint.size(), 13u);
        EXPECT_EQ(dict_res.chiral_restraint.size(), 0u);
        EXPECT_EQ(dict_res.plane_restraint.size(), 0u);
    }

    // Case 3: null pointer
    {
        coot::dictionary_residue_restraints_t dict_res =
            coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t_gemmi(nullptr);
        EXPECT_EQ(dict_res.residue_info.comp_id, "");
        EXPECT_EQ(dict_res.atom_info.size(), 0u);
        EXPECT_EQ(dict_res.bond_restraint.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}