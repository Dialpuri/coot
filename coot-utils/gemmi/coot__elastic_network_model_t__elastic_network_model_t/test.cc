#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <cstdlib>
#include <vector>

TEST(OracleTest, elastic_network_model) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: Select all CA atoms from chain A, build elastic network
    {
        std::vector<const gemmi::Atom*> selected_atoms;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name != "A")
                    continue;
                for (gemmi::Residue& res : chain.residues) {
                    for (gemmi::Atom& atom : res.atoms) {
                        if (atom.name == "CA") {
                            selected_atoms.push_back(&atom);
                        }
                    }
                }
            }
        }

        EXPECT_EQ(selected_atoms.size(), 267u);

        std::srand(1000);
        coot::elastic_network_model_t_gemmi enm(selected_atoms, 5.0, 15.0, 1000);

        std::vector<coot::elastic_network_item_t_gemmi> &items = enm.d;
        EXPECT_EQ(items.size(), 1000u);

        if (!items.empty()) {
            EXPECT_NEAR(items[0].spring_constant, 0.1, 1e-4);
            EXPECT_EQ(items[0].at_1 ? items[0].at_1->name : "null", std::string("CA"));
            EXPECT_EQ(items[0].at_2 ? items[0].at_2->name : "null", std::string("CA"));
        }
    }

    // Case 2: Very tight distance — should produce no springs
    {
        std::vector<const gemmi::Atom*> selected_atoms;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name != "A")
                    continue;
                for (gemmi::Residue& res : chain.residues) {
                    for (gemmi::Atom& atom : res.atoms) {
                        if (atom.name == "CA") {
                            selected_atoms.push_back(&atom);
                        }
                    }
                }
            }
        }

        EXPECT_EQ(selected_atoms.size(), 267u);

        std::srand(1000);
        coot::elastic_network_model_t_gemmi enm(selected_atoms, 100.0, 200.0, 1000);

        std::vector<coot::elastic_network_item_t_gemmi> &items = enm.d;
        EXPECT_EQ(items.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}