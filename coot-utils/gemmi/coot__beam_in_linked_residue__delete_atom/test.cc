#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, beam_in_linked_residue_delete_atom) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());

    gemmi::Model &model = st.models[0];

    // Get first residue from first chain
    ASSERT_FALSE(model.chains.empty());
    gemmi::Chain &chain = model.chains[0];
    ASSERT_FALSE(chain.residues.empty());

    gemmi::Residue &res = chain.residues[0];

    std::string chain_id = chain.name;
    int seqnum = res.seqid.num.value;
    std::string resname = res.name;

    // Count atoms before
    int n_atoms_before = static_cast<int>(res.atoms.size());

    // Find first atom name
    std::string first_atom_name;
    if (n_atoms_before > 0) {
        first_atom_name = res.atoms[0].name;
    }

    EXPECT_EQ(chain_id, "A");
    EXPECT_EQ(seqnum, 1);
    EXPECT_EQ(resname, "MET");

    // Case 1: delete the first atom
    {
        // Call delete_atom
        EXPECT_NO_THROW(coot::beam_in_linked_residue::delete_atom_gemmi(&res, first_atom_name));

        // Count atoms after
        int n_atoms_after = static_cast<int>(res.atoms.size());

        EXPECT_EQ(n_atoms_before, 8);
        EXPECT_EQ(n_atoms_after, 7);
    }

    // Case 2: try to delete an atom that doesn't exist
    {
        int n_atoms_before2 = static_cast<int>(res.atoms.size());

        std::string fake_atom_name = "ZZZZ";

        EXPECT_EQ(n_atoms_before2, 7);

        EXPECT_NO_THROW(coot::beam_in_linked_residue::delete_atom_gemmi(&res, fake_atom_name));

        int n_atoms_after2 = static_cast<int>(res.atoms.size());

        EXPECT_EQ(n_atoms_after2, 7);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}