#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test: delete hydrogens from ILE 10 (0-indexed residue index 9)
TEST(OracleTest, DeleteHydrogensFromILE10) {
  // Load the PDB file with hydrogens
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

  // case 1: delete hydrogens from residue ILE 10 (index 9 in 0-based)
  {
    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];
    gemmi::Residue& res = chain.residues[9];

    int n_atoms_before = res.atoms.size();
    std::cout << "INPUT residue: " << chain.name << " " << res.seqid.num.value << " " << res.name << std::endl;
    std::cout << "INPUT n_atoms_before: " << n_atoms_before << std::endl;
    EXPECT_EQ(n_atoms_before, 19);

    coot::reduce::delete_atom_by_name_gemmi("", res);

    int n_atoms_after = res.atoms.size();
    int atoms_deleted = n_atoms_before - n_atoms_after;
    std::cout << "OUTPUT n_atoms_after: " << n_atoms_after << std::endl;
    std::cout << "OUTPUT atoms_deleted: " << atoms_deleted << std::endl;
    EXPECT_EQ(n_atoms_after, 8);
    EXPECT_EQ(atoms_deleted, 11);
  }

  // case 2: try again on the same residue (should delete nothing since hydrogens are already gone)
  {
    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];
    gemmi::Residue& res = chain.residues[9];

    int n_atoms_before = res.atoms.size();
    std::cout << "INPUT residue: " << chain.name << " " << res.seqid.num.value << " " << res.name << std::endl;
    std::cout << "INPUT n_atoms_before: " << n_atoms_before << std::endl;
    // Should still be 8 (non-hydrogens only)
    EXPECT_EQ(n_atoms_before, 8);

    coot::reduce::delete_atom_by_name_gemmi("", res);

    int n_atoms_after = res.atoms.size();
    int atoms_deleted = n_atoms_before - n_atoms_after;
    std::cout << "OUTPUT n_atoms_after: " << n_atoms_after << std::endl;
    std::cout << "OUTPUT atoms_deleted: " << atoms_deleted << std::endl;
    // Should still be 8, no additional atoms deleted
    EXPECT_EQ(n_atoms_after, 8);
    EXPECT_EQ(atoms_deleted, 0);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}