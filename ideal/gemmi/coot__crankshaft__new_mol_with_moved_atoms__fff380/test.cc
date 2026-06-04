#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a gemmi::Residue by chain name and sequence number
static gemmi::Residue* get_residue_gemmi(
    gemmi::Model& model,
    const std::string& chain_id,
    int seqnum)
{
    for (auto& chain : model.chains) {
        if (chain.name == chain_id) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == seqnum) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, crankshaft) {
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& mol = st.models[0];

        // Get residues 9-14 from chain A
        std::vector<gemmi::Residue*> residues;
        for (int seq = 9; seq <= 14; seq++) {
            gemmi::Residue* r = get_residue_gemmi(mol, "A", seq);
            if (r) {
                residues.push_back(r);
            }
        }

        EXPECT_EQ(residues.size(), 6u);
        EXPECT_STREQ(residues[0]->name.c_str(), "LYS");
        EXPECT_EQ(residues[0]->seqid.num.value, 9);
        EXPECT_EQ(residues[0]->atoms.size(), 9u);

        EXPECT_STREQ(residues[1]->name.c_str(), "ILE");
        EXPECT_EQ(residues[1]->seqid.num.value, 10);
        EXPECT_EQ(residues[1]->atoms.size(), 8u);

        EXPECT_STREQ(residues[2]->name.c_str(), "GLY");
        EXPECT_EQ(residues[2]->seqid.num.value, 11);
        EXPECT_EQ(residues[2]->atoms.size(), 4u);

        EXPECT_STREQ(residues[3]->name.c_str(), "GLU");
        EXPECT_EQ(residues[3]->seqid.num.value, 12);
        EXPECT_EQ(residues[3]->atoms.size(), 9u);

        EXPECT_STREQ(residues[4]->name.c_str(), "GLY");
        EXPECT_EQ(residues[4]->seqid.num.value, 13);
        EXPECT_EQ(residues[4]->atoms.size(), 4u);

        EXPECT_STREQ(residues[5]->name.c_str(), "THR");
        EXPECT_EQ(residues[5]->seqid.num.value, 14);
        EXPECT_EQ(residues[5]->atoms.size(), 7u);

        // Build crankshaft_set from first 4 residues
        coot::crankshaft_set cs0(*residues[0], *residues[1], *residues[2], *residues[3]);
        EXPECT_EQ(cs0.v.size(), 8u);

        // Check atoms in residue 9 (LYS)
        gemmi::Residue *res9 = residues[0];
        EXPECT_EQ(res9->atoms.size(), 9u);
        std::vector<std::string> atom_names;
        for (std::size_t ai = 0; ai < res9->atoms.size(); ai++) {
            gemmi::Atom& a = res9->atoms[ai];
            // gemmi atom.name is unpadded; produce PDB-style 4-char padded name
            // (leading space + name, right-padded to 4)
            std::string padded = " " + a.name;
            padded.resize(4, ' ');
            atom_names.push_back(padded);
        }
        EXPECT_EQ(atom_names.size(), 9u);
        EXPECT_STREQ(atom_names[0].c_str(), " N  ");
        EXPECT_STREQ(atom_names[1].c_str(), " CA ");
        EXPECT_STREQ(atom_names[2].c_str(), " C  ");
        EXPECT_STREQ(atom_names[3].c_str(), " O  ");
        EXPECT_STREQ(atom_names[4].c_str(), " CB ");
        EXPECT_STREQ(atom_names[5].c_str(), " CG ");
        EXPECT_STREQ(atom_names[6].c_str(), " CD ");
        EXPECT_STREQ(atom_names[7].c_str(), " CE ");
        EXPECT_STREQ(atom_names[8].c_str(), " NZ ");

        // Additional coverage: test ca_1 and ca_2 are set
        EXPECT_NE(cs0.ca_1, nullptr);
        EXPECT_NE(cs0.ca_2, nullptr);

        // Test new_mol_with_moved_atoms_gemmi with triple crankshaft set
        coot::crankshaft::scored_triple_angle_set_t_gemmi sas(
            *residues[0], *residues[1], *residues[2],
            *residues[3], *residues[4], *residues[5],
            {0.1f, 0.2f, 0.3f}
        );

        gemmi::Model mol_new = coot::crankshaft::new_mol_with_moved_atoms_gemmi(mol, sas);

        // Verify the copy has the same structure
        EXPECT_EQ(mol_new.chains.size(), mol.chains.size());

        // Verify original model was restored (atoms in mol should be same as before)
        // Count atoms in original model
        size_t atom_count_orig = 0;
        for (auto& chain : mol.chains) {
            for (auto& res : chain.residues) {
                atom_count_orig += res.atoms.size();
            }
        }
        size_t atom_count_new = 0;
        for (auto& chain : mol_new.chains) {
            for (auto& res : chain.residues) {
                atom_count_new += res.atoms.size();
            }
        }
        EXPECT_EQ(atom_count_orig, atom_count_new);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}