#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, water_coordination_t_add_contacts) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // --- Case 1: water + nearby protein atoms (should find contacts with large distance) ---
    {
        coot::util::water_coordination_t wc;

        // Build water selection: all atoms from HOH residues
        std::vector<gemmi::Atom const*> water_sel;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.name == "HOH") {
                        for (gemmi::Atom& atom : res.atoms) {
                            water_sel.push_back(&atom);
                        }
                    }
                }
            }
        }
        int n_water = static_cast<int>(water_sel.size());

        // Build protein selection: all atoms in chain "A"
        std::vector<gemmi::Atom const*> protein_sel;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue& res : chain.residues) {
                        for (gemmi::Atom& atom : res.atoms) {
                            protein_sel.push_back(&atom);
                        }
                    }
                }
            }
        }
        int n_protein = static_cast<int>(protein_sel.size());

        double my_mat[4][4];
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                my_mat[i][j] = (i == j) ? 1.0 : 0.0;

        unsigned int size_before = wc.atom_contacts.size();

        double min_dist = 0.5;
        double max_dist = 50.0;

        wc.add_contacts_gemmi(water_sel, protein_sel,
                              min_dist, max_dist, my_mat);

        unsigned int size_after = wc.atom_contacts.size();

        EXPECT_EQ(n_water, 1);
        EXPECT_EQ(n_protein, 2106);
        EXPECT_NEAR(min_dist, 0.5, 1e-4);
        EXPECT_NEAR(max_dist, 50.0, 1e-4);
        EXPECT_EQ(size_before, 0u);
        EXPECT_EQ(size_after, 1u);

        if (size_after > 0) {
            std::string padded_name = " " + wc.atom_contacts[0].central->name;
            while (padded_name.size() < 4) padded_name += " ";
            EXPECT_STREQ(padded_name.c_str(), " O  ");
            EXPECT_EQ(wc.atom_contacts[0].contact_atoms.size(), 2106u);
        }
    }

    // --- Case 2: tiny distance (should find no contacts) ---
    {
        coot::util::water_coordination_t wc;

        // Build water selection: all atoms from HOH residues
        std::vector<gemmi::Atom const*> water_sel;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.name == "HOH") {
                        for (gemmi::Atom& atom : res.atoms) {
                            water_sel.push_back(&atom);
                        }
                    }
                }
            }
        }
        int n_water = static_cast<int>(water_sel.size());

        // Build protein selection: all atoms in chain "A"
        std::vector<gemmi::Atom const*> protein_sel;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue& res : chain.residues) {
                        for (gemmi::Atom& atom : res.atoms) {
                            protein_sel.push_back(&atom);
                        }
                    }
                }
            }
        }
        int n_protein = static_cast<int>(protein_sel.size());

        double my_mat[4][4];
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                my_mat[i][j] = (i == j) ? 1.0 : 0.0;

        unsigned int size_before = wc.atom_contacts.size();

        double min_dist = 0.5;
        double max_dist = 0.6;

        wc.add_contacts_gemmi(water_sel, protein_sel,
                              min_dist, max_dist, my_mat);

        unsigned int size_after = wc.atom_contacts.size();

        EXPECT_EQ(n_water, 1);
        EXPECT_EQ(n_protein, 2106);
        EXPECT_NEAR(min_dist, 0.5, 1e-4);
        EXPECT_NEAR(max_dist, 0.6, 1e-4);
        EXPECT_EQ(size_before, 0u);
        EXPECT_EQ(size_after, 0u);
    }

    // --- Case 3: moderate distance (hydrogen bond range) — should find some contacts ---
    {
        coot::util::water_coordination_t wc;

        // Build water selection: all atoms from HOH residues
        std::vector<gemmi::Atom const*> water_sel;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.name == "HOH") {
                        for (gemmi::Atom& atom : res.atoms) {
                            water_sel.push_back(&atom);
                        }
                    }
                }
            }
        }
        int n_water = static_cast<int>(water_sel.size());

        // Build protein selection: all atoms in chain "A"
        std::vector<gemmi::Atom const*> protein_sel;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue& res : chain.residues) {
                        for (gemmi::Atom& atom : res.atoms) {
                            protein_sel.push_back(&atom);
                        }
                    }
                }
            }
        }
        int n_protein = static_cast<int>(protein_sel.size());

        double my_mat[4][4];
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                my_mat[i][j] = (i == j) ? 1.0 : 0.0;

        unsigned int size_before = wc.atom_contacts.size();

        double min_dist = 0.5;
        double max_dist = 3.5;

        wc.add_contacts_gemmi(water_sel, protein_sel,
                              min_dist, max_dist, my_mat);

        unsigned int size_after = wc.atom_contacts.size();

        EXPECT_EQ(n_water, 1);
        EXPECT_EQ(n_protein, 2106);
        EXPECT_NEAR(min_dist, 0.5, 1e-4);
        EXPECT_NEAR(max_dist, 3.5, 1e-4);
        EXPECT_EQ(size_before, 0u);
        EXPECT_GT(size_after, 0u);
        EXPECT_LT(size_after, 2106u);

        if (size_after > 0) {
            std::string padded_name = " " + wc.atom_contacts[0].central->name;
            while (padded_name.size() < 4) padded_name += " ";
            EXPECT_STREQ(padded_name.c_str(), " O  ");
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}