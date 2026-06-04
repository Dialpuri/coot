#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <map>
#include <list>
#include <vector>
#include <utility>
#include <algorithm>
#include "geometry/protein-geometry.hh"
#include "ideal/simple-restraint.hh"

static const gemmi::Residue* find_residue(const gemmi::Model& model, const std::string& chain_id, int resno) {
  for (auto& chain : model.chains) {
    if (chain.name == chain_id) {
      for (auto& res : chain.residues) {
        if (res.seqid.num.value == resno) {
          return &res;
        }
      }
    }
  }
  return nullptr;
}

TEST(OracleTest, is_in_same_ring) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  EXPECT_FALSE(st.models.empty());
  gemmi::Model& model = st.models[0];

  coot::protein_geometry geom;
  geom.init_standard();

  const gemmi::Residue* his_res = find_residue(model, "A", 60);
  EXPECT_NE(his_res, nullptr);
  EXPECT_EQ(his_res->name, "HIS");

  // Case 1: HIS residue - CG and CD2 are in the same ring
  {
    std::map<std::string, std::pair<bool, std::vector<std::list<std::string>>>> cache;
    std::string atom1 = " CG ";
    std::string atom2 = " CD2";
    bool result = coot::is_in_same_ring_gemmi(
        0, *his_res, cache, atom1, atom2, geom);
    EXPECT_TRUE(result);
  }

  // Case 2: HIS residue - CG and CB are NOT in the same ring
  {
    std::map<std::string, std::pair<bool, std::vector<std::list<std::string>>>> cache;
    std::string atom1 = " CG ";
    std::string atom2 = " CB ";
    bool result = coot::is_in_same_ring_gemmi(
        0, *his_res, cache, atom1, atom2, geom);
    EXPECT_FALSE(result);
  }

  // Case 3: MET residue - CA and CB (no ring)
  {
    const gemmi::Residue* met_res = find_residue(model, "A", 1);
    EXPECT_NE(met_res, nullptr);
    EXPECT_EQ(met_res->name, "MET");
    std::map<std::string, std::pair<bool, std::vector<std::list<std::string>>>> cache;
    std::string atom1 = " CA ";
    std::string atom2 = " CB ";
    bool result = coot::is_in_same_ring_gemmi(
        0, *met_res, cache, atom1, atom2, geom);
    EXPECT_FALSE(result);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}