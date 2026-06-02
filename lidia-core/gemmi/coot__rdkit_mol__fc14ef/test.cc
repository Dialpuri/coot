#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"
#include "api/molecules-container.hh"
#include "api/coot-molecule.hh"
#include <stdexcept>

TEST(OracleTest, rdkit_mol_gemmi) {
    molecules_container_t mc;
    mc.geometry_init_standard();

    int imol = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_GE(imol, 0);

    int imol_enc = mc.get_imol_enc_any();
    coot::protein_geometry &pg = mc.get_geometry();

    // Load with gemmi to get CRA
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model &model = st.models[0];

    // Find chain A, residue 21 (ALA)
    gemmi::CRA alala_cra{nullptr, nullptr, nullptr};
    for (auto &chain : model.chains) {
        if (chain.name == "A") {
            for (auto &res : chain.residues) {
                if (res.seqid.num.value == 21 && res.name == "ALA") {
                    alala_cra = gemmi::CRA{&chain, &res, nullptr};
                    break;
                }
            }
        }
    }
    ASSERT_NE(alala_cra.residue, nullptr);

    // Case 1: valid standard residue — expects error due to OXT atom
    {
        std::string res_name = alala_cra.residue->name;
        EXPECT_EQ(res_name, "ALA");

        auto p = pg.get_monomer_restraints_at_least_minimal(res_name, imol_enc);
        ASSERT_TRUE(p.first);
        EXPECT_EQ(static_cast<int>(p.second.atom_info.size()), 13);
        EXPECT_EQ(static_cast<int>(p.second.bond_restraint.size()), 12);

        EXPECT_THROW({
            RDKit::RWMol rdkit_mol = coot::rdkit_mol_gemmi(alala_cra, p.second, "", false);
        }, std::runtime_error);
    }

    // Case 2: null residue — should throw
    {
        gemmi::CRA null_cra{nullptr, nullptr, nullptr};
        coot::dictionary_residue_restraints_t empty_restraints;
        EXPECT_THROW({
            RDKit::RWMol rdkit_mol = coot::rdkit_mol_gemmi(null_cra, empty_restraints, "", false);
        }, std::runtime_error);
    }

    // Case 3: residue not in dictionary — construct fake restraints for unknown residue
    {
        coot::dictionary_residue_restraints_t fake_restraints;
        fake_restraints.residue_info.name = "XXX";

        EXPECT_EQ(static_cast<int>(fake_restraints.atom_info.size()), 0);

        RDKit::RWMol rdkit_mol = coot::rdkit_mol_gemmi(alala_cra, fake_restraints, "", false);
        EXPECT_EQ(rdkit_mol.getNumAtoms(), 0u);
        EXPECT_EQ(rdkit_mol.getNumBonds(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}