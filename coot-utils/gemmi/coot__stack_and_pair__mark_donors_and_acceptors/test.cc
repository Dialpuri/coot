#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, stack_and_pair_mark_donors_and_acceptors) {
    // Load PDB into gemmi structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Create protein geometry with standard dictionary
    coot::protein_geometry geom;
    geom.init_standard();

    int imol_enc = -999999; // IMOL_ENC_ANY

    // --- Case 1: Valid selection (all atoms in chain A) ---
    {
        // Collect CRAs for chain A in model 0
        std::vector<gemmi::CRA> atoms;
        gemmi::Model& model = st.models[0];
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    atoms.push_back({&chain, &res, &atom});
                }
            }
        }

        EXPECT_EQ(static_cast<int>(atoms.size()), 2106);
        EXPECT_EQ(imol_enc, -999999);

        std::vector<int> hb_types = coot::stack_and_pair::mark_donors_and_acceptors_gemmi(atoms, geom, imol_enc);

        // Verify: read back hb_types for first 100 atoms
        int sample_count = 0;
        int donor_count = 0;
        int acceptor_count = 0;
        int both_count = 0;
        int unassigned_count = 0;
        for (size_t i = 0; i < atoms.size() && i < 100; i++) {
            int hb_type = hb_types[i];
            sample_count++;
            if (hb_type == 1) donor_count++;        // HB_DONOR
            else if (hb_type == 2) acceptor_count++; // HB_ACCEPTOR
            else if (hb_type == 3) both_count++;     // HB_BOTH
            else unassigned_count++;                 // HB_UNASSIGNED or other
        }
        EXPECT_EQ(sample_count, 100);
        EXPECT_EQ(donor_count, 16);
        EXPECT_EQ(acceptor_count, 21);
        EXPECT_EQ(both_count, 0);
        EXPECT_EQ(unassigned_count, 63);
    }

    // --- Case 2: Empty selection (no atoms selected - chain Z doesn't exist) ---
    {
        std::vector<gemmi::CRA> atoms;
        gemmi::Model& model = st.models[0];
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "Z") continue;
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    atoms.push_back({&chain, &res, &atom});
                }
            }
        }

        EXPECT_EQ(static_cast<int>(atoms.size()), 0);
        EXPECT_EQ(imol_enc, -999999);

        std::vector<int> hb_types = coot::stack_and_pair::mark_donors_and_acceptors_gemmi(atoms, geom, imol_enc);

        EXPECT_EQ(static_cast<int>(hb_types.size()), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}