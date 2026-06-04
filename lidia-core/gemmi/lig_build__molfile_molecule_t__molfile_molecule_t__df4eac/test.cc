#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <sstream>
#include "function.hh"

// Helper: find residue by chain ID and sequence number in gemmi, returning pointer
static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seq_num) {
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, molfile_molecule_t) {
    {
        // Case 1: With bond restraints
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_EQ(st.models.size(), 1u);
        gemmi::Model& model = st.models[0];

        gemmi::Residue *res = find_residue(model, "A", 1);
        ASSERT_NE(res, nullptr);

        EXPECT_EQ(res->seqid.num.value, 1);
        EXPECT_EQ(res->name, "MET");
        EXPECT_EQ(res->atoms.size(), 8u);

        coot::dictionary_residue_restraints_t restraints;
        restraints.residue_info.comp_id = res->name;

        coot::dict_bond_restraint_t bond1(" N  ", " CA ", "bond_1_2", coot::dict_bond_restraint_t::UNASSIGNED);
        restraints.bond_restraint.push_back(bond1);

        coot::dict_bond_restraint_t bond2(" CA ", " C  ", "bond_1_2", coot::dict_bond_restraint_t::UNASSIGNED);
        restraints.bond_restraint.push_back(bond2);

        EXPECT_EQ(restraints.bond_restraint[0].atom_id_1_4c(), " N  ");
        EXPECT_EQ(restraints.bond_restraint[0].atom_id_2_4c(), " CA ");
        EXPECT_EQ(restraints.bond_restraint[1].atom_id_1_4c(), " CA ");
        EXPECT_EQ(restraints.bond_restraint[1].atom_id_2_4c(), " C  ");
        EXPECT_EQ(restraints.bond_restraint.size(), 2u);

        lig_build::molfile_molecule_t mol = lig_build::molfile_molecule_t_gemmi(res, restraints);

        EXPECT_EQ(mol.atoms.size(), 9u);
        EXPECT_EQ(mol.bonds.size(), 2u);

        // atom[0]: null/blank atom (name is empty, xyz is 0,0,0)
        EXPECT_TRUE(mol.atoms[0].name.empty());
        EXPECT_TRUE(mol.atoms[0].element.empty());
        EXPECT_FALSE(mol.atoms[0].aromatic);
        EXPECT_EQ(mol.atoms[0].chiral, 0);
        EXPECT_EQ(mol.atoms[0].formal_charge, 0);
        EXPECT_NEAR(mol.atoms[0].atom_position.x(), 0.0, 1e-4);
        EXPECT_NEAR(mol.atoms[0].atom_position.y(), 0.0, 1e-4);
        EXPECT_NEAR(mol.atoms[0].atom_position.z(), 0.0, 1e-4);

        // atom[1]: N (name/element are space-padded in raw fields, compare formatted)
        {
            std::stringstream ss; ss << mol.atoms[1];
            std::string s = ss.str();
            EXPECT_NE(s.find("N"), std::string::npos);
            EXPECT_NE(s.find("N:"), std::string::npos);
            EXPECT_EQ(mol.atoms[1].chiral, 0);
            EXPECT_EQ(mol.atoms[1].formal_charge, 0);
            EXPECT_NEAR(mol.atoms[1].atom_position.x(), 12.33, 1e-1);
            EXPECT_NEAR(mol.atoms[1].atom_position.y(), 2.772, 1e-1);
            EXPECT_NEAR(mol.atoms[1].atom_position.z(), 74.95, 1e-1);
        }

        // atom[2]: CA
        {
            std::stringstream ss; ss << mol.atoms[2];
            std::string s = ss.str();
            EXPECT_NE(s.find("CA"), std::string::npos);
            EXPECT_NEAR(mol.atoms[2].atom_position.x(), 13.4, 1e-1);
            EXPECT_NEAR(mol.atoms[2].atom_position.y(), 3.78, 1e-1);
            EXPECT_NEAR(mol.atoms[2].atom_position.z(), 74.7, 1e-1);
        }

        // atom[3]: C
        {
            std::stringstream ss; ss << mol.atoms[3];
            EXPECT_NE(ss.str().find("C"), std::string::npos);
        }

        // atom[4]: O
        {
            std::stringstream ss; ss << mol.atoms[4];
            EXPECT_NE(ss.str().find("O:"), std::string::npos);
        }

        // atom[5]: CB
        {
            std::stringstream ss; ss << mol.atoms[5];
            EXPECT_NE(ss.str().find("CB"), std::string::npos);
        }

        // atom[6]: CG
        {
            std::stringstream ss; ss << mol.atoms[6];
            EXPECT_NE(ss.str().find("CG"), std::string::npos);
        }

        // atom[7]: SD
        {
            std::stringstream ss; ss << mol.atoms[7];
            EXPECT_NE(ss.str().find("SD"), std::string::npos);
            EXPECT_NE(ss.str().find("S:"), std::string::npos);
        }

        // atom[8]: CE
        {
            std::stringstream ss; ss << mol.atoms[8];
            EXPECT_NE(ss.str().find("CE"), std::string::npos);
        }

        // bonds: "1 to 2 type 100" and "2 to 3 type 100"
        EXPECT_EQ(mol.bonds[0].index_1, 1);
        EXPECT_EQ(mol.bonds[0].index_2, 2);
        EXPECT_EQ(static_cast<int>(mol.bonds[0].bond_type), 100);
        EXPECT_EQ(mol.bonds[1].index_1, 2);
        EXPECT_EQ(mol.bonds[1].index_2, 3);
        EXPECT_EQ(static_cast<int>(mol.bonds[1].bond_type), 100);
    }

    {
        // Case 2: Empty restraints — atoms still created, no bonds
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_EQ(st.models.size(), 1u);
        gemmi::Model& model = st.models[0];

        gemmi::Residue *res = find_residue(model, "A", 1);
        ASSERT_NE(res, nullptr);

        coot::dictionary_residue_restraints_t empty_restraints;
        empty_restraints.residue_info.comp_id = res->name;

        lig_build::molfile_molecule_t mol_empty = lig_build::molfile_molecule_t_gemmi(res, empty_restraints);
        EXPECT_EQ(mol_empty.atoms.size(), 9u);
        EXPECT_EQ(mol_empty.bonds.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}