#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"

#include <string>
#include <vector>

namespace {

struct WaterInfo {
  std::string chain;
  std::string rname;
  int seqnum;
  double x, y, z;
};

std::vector<WaterInfo> collect_waters(gemmi::Structure& st) {
  std::vector<WaterInfo> result;
  gemmi::Model& model = st.models[0];
  for (gemmi::Chain& chain : model.chains) {
    for (gemmi::Residue& residue : chain.residues) {
      std::string rname(residue.name);
      if (rname == "WAT" || rname == "HOH") {
        for (gemmi::Atom& atom : residue.atoms) {
          WaterInfo wi;
          wi.chain = chain.name;
          wi.rname = rname;
          wi.seqnum = residue.seqid.num.value;
          wi.x = atom.pos.x;
          wi.y = atom.pos.y;
          wi.z = atom.pos.z;
          result.push_back(wi);
        }
      }
    }
  }
  return result;
}

} // namespace

TEST(OracleTest, move_waters_around_protein) {
  // Case 1: Normal case — load PDB, add waters, call function
  {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_GE(st.models.size(), 1);

    // Set cell and spacegroup (required for symmetry_move_atoms to work)
    st.cell = gemmi::UnitCell(40.0, 40.0, 40.0, 90.0, 90.0, 90.0);
    st.spacegroup_hm = "P 1";

    // Add water residues to the first chain so we can test the function
    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];

    // Create HOH residue at seqnum 900
    gemmi::Residue hoh1;
    hoh1.name = "HOH";
    hoh1.seqid.num.value = 900;
    hoh1.seqid.icode = ' ';
    hoh1.het_flag = 'H';
    gemmi::Atom o1;
    o1.name = "O";
    o1.element = gemmi::Element("O");
    o1.pos = gemmi::Position(1.0, 1.0, 1.0);
    o1.occ = 1.0;
    o1.b_iso = 0.0;
    hoh1.atoms.push_back(o1);
    chain.residues.push_back(hoh1);

    // Create HOH residue at seqnum 901
    gemmi::Residue hoh2;
    hoh2.name = "HOH";
    hoh2.seqid.num.value = 901;
    hoh2.seqid.icode = ' ';
    hoh2.het_flag = 'H';
    gemmi::Atom o2;
    o2.name = "O";
    o2.element = gemmi::Element("O");
    o2.pos = gemmi::Position(2.0, 2.0, 2.0);
    o2.occ = 1.0;
    o2.b_iso = 0.0;
    hoh2.atoms.push_back(o2);
    chain.residues.push_back(hoh2);

    // Create WAT residue at seqnum 902
    gemmi::Residue wat;
    wat.name = "WAT";
    wat.seqid.num.value = 902;
    wat.seqid.icode = ' ';
    wat.het_flag = 'H';
    gemmi::Atom o3;
    o3.name = "O";
    o3.element = gemmi::Element("O");
    o3.pos = gemmi::Position(3.0, 3.0, 3.0);
    o3.occ = 1.0;
    o3.b_iso = 0.0;
    wat.atoms.push_back(o3);
    chain.residues.push_back(wat);

    // Collect water positions before
    std::vector<WaterInfo> waters_before = collect_waters(st);
    EXPECT_EQ(waters_before.size(), 4u);

    // Verify initial water positions (oracle ground truth)
    EXPECT_EQ(waters_before[0].chain, "A");
    EXPECT_EQ(waters_before[0].rname, "HOH");
    EXPECT_EQ(waters_before[0].seqnum, 900);
    EXPECT_NEAR(waters_before[0].x, 1.0, 1e-4);
    EXPECT_NEAR(waters_before[0].y, 1.0, 1e-4);
    EXPECT_NEAR(waters_before[0].z, 1.0, 1e-4);

    EXPECT_EQ(waters_before[1].chain, "A");
    EXPECT_EQ(waters_before[1].rname, "HOH");
    EXPECT_EQ(waters_before[1].seqnum, 901);
    EXPECT_NEAR(waters_before[1].x, 2.0, 1e-4);
    EXPECT_NEAR(waters_before[1].y, 2.0, 1e-4);
    EXPECT_NEAR(waters_before[1].z, 2.0, 1e-4);

    EXPECT_EQ(waters_before[2].chain, "A");
    EXPECT_EQ(waters_before[2].rname, "WAT");
    EXPECT_EQ(waters_before[2].seqnum, 902);
    EXPECT_NEAR(waters_before[2].x, 3.0, 1e-4);
    EXPECT_NEAR(waters_before[2].y, 3.0, 1e-4);
    EXPECT_NEAR(waters_before[2].z, 3.0, 1e-4);

    EXPECT_EQ(waters_before[3].chain, "B");
    EXPECT_EQ(waters_before[3].rname, "HOH");
    EXPECT_EQ(waters_before[3].seqnum, 1);
    EXPECT_NEAR(waters_before[3].x, 30.744, 1e-3);
    EXPECT_NEAR(waters_before[3].y, 8.475, 1e-3);
    EXPECT_NEAR(waters_before[3].z, 36.748, 1e-3);

    int n_moved = coot::util::move_waters_around_protein_gemmi(st);
    EXPECT_EQ(n_moved, 4);

    // Collect water positions after
    std::vector<WaterInfo> waters_after = collect_waters(st);
    ASSERT_EQ(waters_after.size(), 4u);

    // Verify final water positions (oracle ground truth)
    EXPECT_EQ(waters_after[0].chain, "A");
    EXPECT_EQ(waters_after[0].rname, "HOH");
    EXPECT_EQ(waters_after[0].seqnum, 900);
    EXPECT_NEAR(waters_after[0].x, 41.0, 1e-4);
    EXPECT_NEAR(waters_after[0].y, 1.0, 1e-4);
    EXPECT_NEAR(waters_after[0].z, 41.0, 1e-4);

    EXPECT_EQ(waters_after[1].chain, "A");
    EXPECT_EQ(waters_after[1].rname, "HOH");
    EXPECT_EQ(waters_after[1].seqnum, 901);
    EXPECT_NEAR(waters_after[1].x, 42.0, 1e-4);
    EXPECT_NEAR(waters_after[1].y, 2.0, 1e-4);
    EXPECT_NEAR(waters_after[1].z, 42.0, 1e-4);

    EXPECT_EQ(waters_after[2].chain, "A");
    EXPECT_EQ(waters_after[2].rname, "WAT");
    EXPECT_EQ(waters_after[2].seqnum, 902);
    EXPECT_NEAR(waters_after[2].x, 43.0, 1e-4);
    EXPECT_NEAR(waters_after[2].y, 3.0, 1e-4);
    EXPECT_NEAR(waters_after[2].z, 43.0, 1e-4);

    // Water in chain B should stay in place (already far from protein)
    EXPECT_EQ(waters_after[3].chain, "B");
    EXPECT_EQ(waters_after[3].rname, "HOH");
    EXPECT_EQ(waters_after[3].seqnum, 1);
    EXPECT_NEAR(waters_after[3].x, 30.744, 1e-3);
    EXPECT_NEAR(waters_after[3].y, 8.475, 1e-3);
    EXPECT_NEAR(waters_after[3].z, 36.748, 1e-3);
  }

  // Case 2: Load PDB without adding extra waters — tests with only existing waters
  {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(st.models.size(), 1);

    st.cell = gemmi::UnitCell(40.0, 40.0, 40.0, 90.0, 90.0, 90.0);
    st.spacegroup_hm = "P 1";

    // Verify only the existing water from the PDB (chain B, HOH 1)
    std::vector<WaterInfo> waters_before = collect_waters(st);
    EXPECT_EQ(waters_before.size(), 1u);
    EXPECT_EQ(waters_before[0].chain, "B");
    EXPECT_EQ(waters_before[0].rname, "HOH");
    EXPECT_EQ(waters_before[0].seqnum, 1);

    int n_moved = coot::util::move_waters_around_protein_gemmi(st);
    // With one water, the function processes it (returns 1 regardless of actual displacement)
    EXPECT_EQ(n_moved, 1);

    // Verify water position (may or may not change depending on proximity to protein)
    std::vector<WaterInfo> waters_after = collect_waters(st);
    EXPECT_EQ(waters_after.size(), 1u);
    EXPECT_NEAR(waters_after[0].x, 30.744, 1e-3);
    EXPECT_NEAR(waters_after[0].y, 8.475, 1e-3);
    EXPECT_NEAR(waters_after[0].z, 36.748, 1e-3);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}