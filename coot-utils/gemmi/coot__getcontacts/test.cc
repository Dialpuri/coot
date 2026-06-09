#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, getcontacts) {
    coot::protein_geometry geom;
    geom.init_standard();

    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];
    gemmi::UnitCell cell = st.cell;

    // Case 1: residue //A/10 (ILE) with geometry — 4-param variant
    {
        std::vector<gemmi::Atom*> atoms;
        std::string monomer_type;

        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value != 10) continue;
                monomer_type = res.name;
                for (gemmi::Atom& atom : res.atoms) {
                    atoms.push_back(&atom);
                }
                break;
            }
        }

        EXPECT_FALSE(atoms.empty());

        contact_info ci = coot::getcontacts_gemmi(model, cell, atoms, monomer_type, 0, &geom);
        EXPECT_EQ(ci.n_contacts(), 7);

        std::vector<std::vector<int>> indices = ci.get_contact_indices_with_reverse_contacts();
        EXPECT_EQ(indices.size(), 8u);
    }

    // Case 2: residue //A/50 (ARG) with geometry — 4-param variant
    {
        std::vector<gemmi::Atom*> atoms;
        std::string monomer_type;

        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value != 50) continue;
                monomer_type = res.name;
                for (gemmi::Atom& atom : res.atoms) {
                    atoms.push_back(&atom);
                }
                break;
            }
        }

        EXPECT_FALSE(atoms.empty());

        contact_info ci = coot::getcontacts_gemmi(model, cell, atoms, monomer_type, 0, &geom);
        EXPECT_EQ(ci.n_contacts(), 10);

        std::vector<std::vector<int>> indices = ci.get_contact_indices_with_reverse_contacts();
        EXPECT_EQ(indices.size(), 11u);
    }

    // Case 3: residue with alt_conf parameter — 5-param variant
    // In the gemmi port, the 5-param variant doesn't exist as a separate signature
    // (alt_conf is not used in the gemmi port). So we test the 4-param variant with same data.
    {
        std::vector<gemmi::Atom*> atoms;
        std::string monomer_type;

        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value != 10) continue;
                monomer_type = res.name;
                for (gemmi::Atom& atom : res.atoms) {
                    atoms.push_back(&atom);
                }
                break;
            }
        }

        EXPECT_FALSE(atoms.empty());

        contact_info ci = coot::getcontacts_gemmi(model, cell, atoms, monomer_type, 0, &geom);
        EXPECT_EQ(ci.n_contacts(), 7);

        std::vector<std::vector<int>> indices = ci.get_contact_indices_with_reverse_contacts();
        EXPECT_EQ(indices.size(), 8u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
