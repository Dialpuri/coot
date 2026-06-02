#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

typedef std::pair<gemmi::Atom*, gemmi::Residue*> sel_entry;

static std::vector<sel_entry>
build_residue_selection(gemmi::Structure& st, char chain_id, int seq_num) {
  std::vector<sel_entry> sel;
  for (auto& model : st.models) {
    for (auto& chain : model.chains) {
      if (chain.name.size() == 0 || chain.name[0] != chain_id) continue;
      for (auto& res : chain.residues) {
        if (res.seqid.num.value != seq_num) continue;
        for (auto& atom : res.atoms) {
          sel.push_back({&atom, &res});
        }
      }
    }
  }
  return sel;
}

static std::vector<sel_entry>
build_atom_selection(gemmi::Structure& st, char chain_id, int seq_num,
                     const std::vector<std::string>& atom_names) {
  std::vector<sel_entry> sel;
  for (auto& model : st.models) {
    for (auto& chain : model.chains) {
      if (chain.name.size() == 0 || chain.name[0] != chain_id) continue;
      for (auto& res : chain.residues) {
        if (res.seqid.num.value != seq_num) continue;
        for (auto& atom : res.atoms) {
          for (const auto& target : atom_names) {
            if (atom.name == target) {
              sel.push_back({&atom, &res});
              break;
            }
          }
        }
      }
    }
  }
  return sel;
}

TEST(OracleTest, mergeable_atom_selections) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

  // Case 1: Two selections on the same non-existent residue (oracle: false)
  {
    auto sel1 = build_residue_selection(st, 'A', 10);
    auto sel2 = build_residue_selection(st, 'A', 10);
    auto result = coot::mergeable_atom_selections_gemmi(sel1, sel2);
    EXPECT_FALSE(result.first);
  }

  // Case 2: Two selections on residues where one doesn't exist (oracle: false)
  {
    auto sel3 = build_residue_selection(st, 'A', 1);
    auto sel4 = build_residue_selection(st, 'B', 50);
    auto result2 = coot::mergeable_atom_selections_gemmi(sel3, sel4);
    EXPECT_FALSE(result2.first);
  }

  // Case 3: Two selections on non-existent adjacent residues (oracle: false)
  {
    auto sel5 = build_residue_selection(st, 'A', 10);
    auto sel6 = build_residue_selection(st, 'A', 11);
    auto result3 = coot::mergeable_atom_selections_gemmi(sel5, sel6);
    EXPECT_FALSE(result3.first);
  }

  // Case 4 (NEW): Two selections on the same existing residue via atom specs
  // Self-contacts skipped by at_1==at_2 guard, cross-atom pairs have different names
  // so no matching pairs → 0 pairs → false
  {
    auto sel7 = build_atom_selection(st, 'A', 1, {"N ", "CA", "C ", "CB"});
    auto sel8 = build_atom_selection(st, 'A', 1, {"N ", "CA", "C ", "CB"});
    auto result4 = coot::mergeable_atom_selections_gemmi(sel7, sel8);
    EXPECT_FALSE(result4.first);
  }

  // Case 5 (NEW): Two HOH selections — O atoms excluded by the function
  {
    auto sel9 = build_atom_selection(st, 'B', 1, {"O  "});
    auto sel10 = build_atom_selection(st, 'B', 1, {"O  "});
    auto result5 = coot::mergeable_atom_selections_gemmi(sel9, sel10);
    EXPECT_FALSE(result5.first);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}