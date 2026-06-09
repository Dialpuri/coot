#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

// Helper: find residue by chain name and sequence number in the first model
static gemmi::Residue* find_residue(gemmi::Structure& st, const std::string& chain_name, int seq_num) {
    if (st.models.empty()) return nullptr;
    gemmi::Model& model = st.models[0];
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_name) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, getcontacts) {
    // Load the PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Initialize protein geometry (REQUIRED)
    coot::protein_geometry geom;
    geom.init_standard();

    int imol = 0;

    // --- Case 1: valid residue (ILE) ---
    {
        gemmi::Residue* res1 = find_residue(st, "A", 10);
        ASSERT_NE(res1, nullptr);

        std::vector<gemmi::Atom*> atoms1;
        for (gemmi::Atom& atom : res1->atoms) {
            atoms1.push_back(&atom);
        }

        int n_atoms1 = static_cast<int>(atoms1.size());
        std::string alt_conf1("");
        std::string monomer_type1(res1->name);

        coot::contact_info::contact_info_result_gemmi ci1 = coot::getcontacts_gemmi(atoms1, alt_conf1, monomer_type1, imol, &geom);

        EXPECT_EQ(monomer_type1, "ILE");
        EXPECT_EQ(n_atoms1, 8);
        EXPECT_EQ(ci1.n_contacts(), 7);
    }

    // --- Case 2: different residue (GLY) ---
    {
        gemmi::Residue* res2 = find_residue(st, "A", 11);
        ASSERT_NE(res2, nullptr);

        std::vector<gemmi::Atom*> atoms2;
        for (gemmi::Atom& atom : res2->atoms) {
            atoms2.push_back(&atom);
        }

        int n_atoms2 = static_cast<int>(atoms2.size());
        std::string alt_conf2("");
        std::string monomer_type2(res2->name);

        coot::contact_info::contact_info_result_gemmi ci2 = coot::getcontacts_gemmi(atoms2, alt_conf2, monomer_type2, imol, &geom);

        EXPECT_EQ(monomer_type2, "GLY");
        EXPECT_EQ(n_atoms2, 4);
        EXPECT_EQ(ci2.n_contacts(), 3);
    }

    // --- Case 3: another residue in chain A (GLU) ---
    {
        gemmi::Residue* res3 = find_residue(st, "A", 12);
        ASSERT_NE(res3, nullptr);

        std::vector<gemmi::Atom*> atoms3;
        for (gemmi::Atom& atom : res3->atoms) {
            atoms3.push_back(&atom);
        }

        int n_atoms3 = static_cast<int>(atoms3.size());
        std::string alt_conf3("");
        std::string monomer_type3(res3->name);

        coot::contact_info::contact_info_result_gemmi ci3 = coot::getcontacts_gemmi(atoms3, alt_conf3, monomer_type3, imol, &geom);

        EXPECT_EQ(monomer_type3, "GLU");
        EXPECT_EQ(n_atoms3, 9);
        EXPECT_EQ(ci3.n_contacts(), 8);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
