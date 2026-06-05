#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, my_atom_colour_map_t_fill_chain_id_map) {
    // Case: fill_chain_id_map with two chains (A, B)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        ASSERT_GT(static_cast<int>(st.models.size()), 0);

        const gemmi::Model &model = st.models[0];

        // Count chains before we use the map
        int n_chains = static_cast<int>(model.chains.size());
        EXPECT_EQ(n_chains, 2);

        // Print chain IDs
        std::string chain_ids;
        for (const gemmi::Chain &ch : model.chains) {
            if (!chain_ids.empty()) chain_ids += ", ";
            chain_ids += ch.name;
        }
        EXPECT_EQ(chain_ids, "A, B");

        // Count selected atoms (all atoms in model 0)
        int n_selected_atoms = 0;
        for (const gemmi::Chain &ch : model.chains)
            for (const gemmi::Residue &res : ch.residues)
                n_selected_atoms += static_cast<int>(res.atoms.size());
        EXPECT_EQ(n_selected_atoms, 2107);

        // Create the colour map and call fill_chain_id_map
        coot::my_atom_colour_map_t acm;

        // BEFORE: atom_colour_map starts with 50 entries of "---"
        unsigned int size_before = acm.atom_colour_map.size();
        EXPECT_EQ(size_before, 50u);

        acm.fill_chain_id_map_gemmi(model);

        // AFTER: atom_colour_map should have entries appended for each chain
        unsigned int size_after = acm.atom_colour_map.size();
        EXPECT_EQ(size_after, 52u);

        // Print the appended chain entries (after the initial 50)
        EXPECT_EQ(acm.atom_colour_map[50], "A");
        EXPECT_EQ(acm.atom_colour_map[51], "B");

        // Test index_for_chain for known chains and an unknown chain
        for (const gemmi::Chain &ch : model.chains) {
            std::string cid = ch.name;
            unsigned int idx = static_cast<unsigned int>(acm.index_for_chain(cid));
            if (cid == "A") {
                EXPECT_EQ(idx, 50u);
            } else if (cid == "B") {
                EXPECT_EQ(idx, 51u);
            }
        }

        // Unknown chain
        unsigned int idx_unknown = static_cast<unsigned int>(acm.index_for_chain("Z"));
        EXPECT_EQ(idx_unknown, 52u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}