#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ---- Tests ----

TEST(OracleTest, assignCharge) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];

    // Case 1: Select all atoms in chain A — scan for atoms with known charges
    {
        // Build selection of all atoms in chain A
        std::vector<gemmi::CRA> selection;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& residue : chain.residues) {
                    for (gemmi::Atom& atom : residue.atoms) {
                        selection.push_back(gemmi::CRA{&chain, &residue, &atom});
                    }
                }
            }
        }

        EXPECT_EQ(selection.size(), 2106u);

        // Scan selection for atoms with known non-zero charges in the table
        std::string charged_names[] = {" OT1", " OT2", " OD1", " OD2", " OE1", " OE2", " NH1", " NH2", " NZ "};
        std::vector<gemmi::CRA*> charged_atoms;
        for (size_t i = 0; i < selection.size() && charged_atoms.size() < 4; i++) {
            gemmi::CRA* cra = &selection[i];
            std::string aname = cra->atom->padded_name();
            for (size_t c = 0; c < sizeof(charged_names)/sizeof(charged_names[0]); c++) {
                if (aname == charged_names[c]) {
                    charged_atoms.push_back(cra);
                    break;
                }
            }
        }

        EXPECT_EQ(charged_atoms.size(), 4u);

        // Verify atom identities match oracle
        EXPECT_EQ(charged_atoms[0]->atom->padded_name(), " OE1");
        EXPECT_EQ(charged_atoms[0]->residue->name, "GLU");
        EXPECT_EQ(charged_atoms[0]->residue->seqid.num.value, 2);

        EXPECT_EQ(charged_atoms[1]->atom->padded_name(), " OE2");
        EXPECT_EQ(charged_atoms[1]->residue->name, "GLU");
        EXPECT_EQ(charged_atoms[1]->residue->seqid.num.value, 2);

        EXPECT_EQ(charged_atoms[2]->atom->padded_name(), " OD1");
        EXPECT_EQ(charged_atoms[2]->residue->name, "ASN");
        EXPECT_EQ(charged_atoms[2]->residue->seqid.num.value, 3);

        EXPECT_EQ(charged_atoms[3]->atom->padded_name(), " OE1");
        EXPECT_EQ(charged_atoms[3]->residue->name, "GLN");
        EXPECT_EQ(charged_atoms[3]->residue->seqid.num.value, 5);

        // Call assignCharge
        CXXChargeTable theChargeTable;
        std::map<const gemmi::Atom*, double> charges_out;
        int result = CXXUtils::assignCharge_gemmi(selection, &theChargeTable, charges_out);
        EXPECT_EQ(result, 0);

        // Check AFTER charges from the output map
        EXPECT_NEAR(charges_out[charged_atoms[0]->atom], -0.5, 1e-4);
        EXPECT_NEAR(charges_out[charged_atoms[1]->atom], -0.5, 1e-4);
        EXPECT_NEAR(charges_out[charged_atoms[2]->atom], 0.0, 1e-4);
        EXPECT_NEAR(charges_out[charged_atoms[3]->atom], 0.0, 1e-4);

        // Charge table lookup examples
        EXPECT_NEAR(theChargeTable.getCharge("ALA", " OT1"), -0.5, 1e-4);
        EXPECT_NEAR(theChargeTable.getCharge("ARG", " NH1"), 0.5, 1e-4);
        EXPECT_NEAR(theChargeTable.getCharge("XXX", "XX"), 0.0, 1e-4);
    }

    // Case 2: Empty selection — verifies guarded path (no atoms to modify)
    {
        // Build selection of all atoms in chain Z (doesn't exist)
        std::vector<gemmi::CRA> selection;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "Z") {
                for (gemmi::Residue& residue : chain.residues) {
                    for (gemmi::Atom& atom : residue.atoms) {
                        selection.push_back(gemmi::CRA{&chain, &residue, &atom});
                    }
                }
            }
        }

        EXPECT_EQ(selection.size(), 0u);

        CXXChargeTable theChargeTable;
        std::map<const gemmi::Atom*, double> charges_out;
        int result = CXXUtils::assignCharge_gemmi(selection, &theChargeTable, charges_out);
        EXPECT_EQ(result, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}