#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <coot/utils/coot-utils.hh>
#include "/lmb/home/jdialpuri/Development/coot-dev/coot/ligand/gemmi/coot__ideal_rna__get_standard_residue_instance/function.hh"

// Helper: check if a residue has an O2' or O2* atom
static bool has_o2_prime(const gemmi::Residue *res) {
  if (!res) return false;
  for (const auto &atom : res->atoms) {
    std::string an = atom.name;
    if (an == "O2'" || an == "O2*") {
      return true;
    }
  }
  return false;
}

TEST(OracleTest, delete_o2_prime) {
  // Load standard residues
  gemmi::Structure std_struct = gemmi::read_pdb_file(coot::package_data_dir() + "/standard-residues.pdb");
  const gemmi::Model& std_model = std_struct.models[0];

  // Case 1: RNA G residue (should have O2' atom)
  {
    gemmi::Residue *res = coot::ideal_rna::get_standard_residue_instance_gemmi("G", std_model);
    ASSERT_NE(res, nullptr);

    // Count atoms before
    int natoms_before = static_cast<int>(res->atoms.size());

    // Check if O2' exists before
    bool has_o2_prime_before = has_o2_prime(res);

    EXPECT_STREQ(res->name.c_str(), "G");
    EXPECT_EQ(natoms_before, 23);
    EXPECT_TRUE(has_o2_prime_before);

    // Call delete_o2_prime
    coot::ideal_rna::delete_o2_prime_gemmi(res);

    // Count atoms after
    int natoms_after = static_cast<int>(res->atoms.size());

    // Check if O2' still exists after
    bool has_o2_prime_after = has_o2_prime(res);

    EXPECT_EQ(natoms_after, 22);
    EXPECT_EQ((natoms_before - natoms_after), 1);
    EXPECT_FALSE(has_o2_prime_after);
  }

  // Case 2: DNA DG residue (should NOT have O2' atom already)
  {
    gemmi::Residue *res = coot::ideal_rna::get_standard_residue_instance_gemmi("DG", std_model);
    ASSERT_NE(res, nullptr);

    int natoms_before = static_cast<int>(res->atoms.size());

    bool has_o2_prime_before = has_o2_prime(res);

    EXPECT_STREQ(res->name.c_str(), "DG");
    EXPECT_EQ(natoms_before, 22);
    EXPECT_FALSE(has_o2_prime_before);

    // Call delete_o2_prime (should be a no-op since DNA lacks O2')
    coot::ideal_rna::delete_o2_prime_gemmi(res);

    int natoms_after = static_cast<int>(res->atoms.size());

    EXPECT_EQ(natoms_after, 22);
    EXPECT_EQ((natoms_before - natoms_after), 0);
  }

  // Case 3: NULL residue (edge case — should do nothing)
  {
    gemmi::Residue *res = nullptr;
    EXPECT_NO_THROW(coot::ideal_rna::delete_o2_prime_gemmi(res));
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}