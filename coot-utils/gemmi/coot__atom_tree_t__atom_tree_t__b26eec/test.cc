#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "geometry/protein-geometry.hh"

// Helper: find a residue by chain name and sequence number in a gemmi Structure
static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seq_num) {
  for (auto& chain : model.chains) {
    if (chain.name != chain_id) continue;
    for (auto& res : chain.residues) {
      if (res.seqid.num.value == seq_num) {
        return &res;
      }
    }
  }
  return nullptr;
}

TEST(OracleTest, atom_tree_construction) {
  // Setup: load molecule with gemmi
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  EXPECT_EQ(st.models.size(), 1u);

  gemmi::Model& model = st.models[0];

  coot::protein_geometry geom;
  geom.init_standard();

  // Case 1: Standard residue (ILE) with restraints
  {
    // cid = "//A/10" means chain A, residue 10
    gemmi::Residue* res = find_residue(model, "A", 10);
    EXPECT_NE(res, nullptr);

    std::string res_name(res->name);
    EXPECT_EQ(res_name, "ILE");

    std::string alt_conf("");

    std::pair<bool, coot::dictionary_residue_restraints_t> restraints_info =
        geom.get_monomer_restraints(res_name, 0);

    EXPECT_TRUE(restraints_info.first);

    bool constructed = false;
    std::string error_msg;
    try {
      coot::atom_tree_t tree;
      tree.atom_tree_t_gemmi(restraints_info.second, *res, alt_conf);
      constructed = true;
    } catch (const std::runtime_error &rte) {
      error_msg = rte.what();
    }
    EXPECT_TRUE(constructed);
  }

  // Case 2: Another residue (ARG) with restraints
  {
    // cid = "//A/50" means chain A, residue 50
    gemmi::Residue* res = find_residue(model, "A", 50);
    EXPECT_NE(res, nullptr);

    std::string res_name(res->name);
    EXPECT_EQ(res_name, "ARG");

    std::string alt_conf("");

    std::pair<bool, coot::dictionary_residue_restraints_t> restraints_info =
        geom.get_monomer_restraints(res_name, 0);

    EXPECT_TRUE(restraints_info.first);

    bool constructed = false;
    try {
      coot::atom_tree_t tree;
      tree.atom_tree_t_gemmi(restraints_info.second, *res, alt_conf);
      constructed = true;
    } catch (const std::runtime_error &rte) {
      // expected for some edge cases
    }
    EXPECT_TRUE(constructed);
  }

  // Case 3: Non-existent residue - verifies error path
  {
    // cid = "//Z/9999" means chain Z, residue 9999
    gemmi::Residue* res = find_residue(model, "Z", 9999);
    EXPECT_EQ(res, nullptr);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}