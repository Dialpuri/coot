#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, contact_info_gemmi) {
    // Load the PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Initialize protein geometry
    coot::protein_geometry geom;
    geom.init_standard();

    int imol = 0;

    // --- Case 1: Select a real residue (chain A, res 1) ---
    {
        gemmi::Model& model = st.models[0];

        // Find chain A, residue 1 (MET)
        gemmi::Residue* res = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& r : chain.residues) {
                    if (r.seqid.num.value == 1) {
                        res = &r;
                        break;
                    }
                }
                break;
            }
        }
        ASSERT_NE(res, nullptr);

        std::string monomer_type(res->name);
        int n_atoms = static_cast<int>(res->atoms.size());

        EXPECT_EQ(monomer_type, "MET");
        EXPECT_EQ(n_atoms, 8);

        // Build atom selection from this residue's atoms
        std::vector<gemmi::Atom*> atom_selection;
        for (gemmi::Atom& atom : res->atoms) {
            atom_selection.push_back(&atom);
        }

        coot::contact_info::contact_info_result ci =
            coot::contact_info::contact_info_gemmi(atom_selection, monomer_type, imol, &geom);

        EXPECT_EQ(ci.n_contacts(), 7);
    }

    // --- Case 2: Empty selection (0 atoms) — edge case ---
    {
        std::vector<gemmi::Atom*> atom_selection;

        coot::contact_info::contact_info_result ci =
            coot::contact_info::contact_info_gemmi(atom_selection, "ALA", imol, &geom);

        EXPECT_EQ(ci.n_contacts(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}