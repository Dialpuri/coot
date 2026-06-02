#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <string>
#include <cstring>

static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int res_seqnum) {
  for (auto& chain : model.chains) {
    if (chain.name == chain_id) {
      for (auto& residue : chain.residues) {
        if (residue.seqid.num.value == res_seqnum) {
          return &residue;
        }
      }
    }
  }
  return nullptr;
}

TEST(OracleTest, residue_atoms_segid) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  gemmi::Model& model = st.models[0];

  // Case 1: normal residue with consistent segids
  {
    gemmi::Residue* res = find_residue(model, "A", 10);
    ASSERT_NE(res, nullptr);
    EXPECT_EQ(res->name, "ILE");

    bool threw = false;
    std::string segid;
    try {
      segid = coot::residue_atoms_segid_gemmi(*res);
    } catch (const std::runtime_error &) {
      threw = true;
    }
    EXPECT_FALSE(threw);
    EXPECT_EQ(segid, "");
  }

  // Case 2: another residue to confirm consistency
  {
    gemmi::Residue* res = find_residue(model, "A", 50);
    ASSERT_NE(res, nullptr);
    EXPECT_EQ(res->name, "ARG");

    bool threw = false;
    std::string segid;
    try {
      segid = coot::residue_atoms_segid_gemmi(*res);
    } catch (const std::runtime_error &) {
      threw = true;
    }
    EXPECT_FALSE(threw);
    EXPECT_EQ(segid, "");
  }

  // Case 3: HOH residue from chain B
  {
    gemmi::Residue* res = find_residue(model, "B", 1);
    ASSERT_NE(res, nullptr);
    EXPECT_EQ(res->name, "HOH");

    bool threw = false;
    std::string segid;
    try {
      segid = coot::residue_atoms_segid_gemmi(*res);
    } catch (const std::runtime_error &) {
      threw = true;
    }
    EXPECT_FALSE(threw);
    EXPECT_EQ(segid, "");
  }

  // Case 4: non-empty segid — set segid on residue to verify core logic
  // (in gemmi, segid is per-residue, not per-atom)
  {
    gemmi::Residue* res = find_residue(model, "A", 10);
    ASSERT_NE(res, nullptr);

    // Set residue segid to "AAA"
    res->segment = "AAA";

    std::string segid = coot::residue_atoms_segid_gemmi(*res);
    EXPECT_EQ(segid, "AAA");
  }

  // Case 5: inconsistent segids should throw
  // In gemmi, segid is per-residue so this cannot happen naturally.
  // We cannot set per-atom segids. The test verifies the original
  // expectation that inconsistent segids would throw, but in gemmi's
  // data model this scenario is impossible within a single residue.
  {
    gemmi::Residue* res = find_residue(model, "A", 50);
    ASSERT_NE(res, nullptr);
    ASSERT_GE(res->atoms.size(), 2u);

    // In gemmi, segid is at the Residue level (inherited from ResidueId),
    // so all atoms in a residue share the same segid by design.
    // This means inconsistent segids within a residue cannot occur.
    // We verify the function does NOT throw (consistent segid across all atoms).
    bool threw = false;
    try {
      coot::residue_atoms_segid_gemmi(*res);
    } catch (const std::runtime_error &) {
      threw = true;
    }
    // In gemmi model, segid is per-residue, so consistency is guaranteed.
    EXPECT_FALSE(threw);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}