#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_contact_indices_for_PRO_residue) {
  // Initialize geometry
  coot::protein_geometry geom;
  geom.init_standard();

  // Load PDB with gemmi
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  ASSERT_FALSE(st.models.empty());
  gemmi::Model& model = st.models[0];

  // Case 1: PRO residue at //A/61
  {
    std::string cid = "//A/61";

    // Find PRO residue at chain A, seqnum 61
    gemmi::Residue* res = nullptr;
    for (gemmi::Chain& ch : model.chains) {
      if (ch.name == "A") {
        for (gemmi::Residue& r : ch.residues) {
          if (r.seqid.num.value == 61 && r.name == "PRO") {
            res = &r;
            break;
          }
        }
      }
    }
    ASSERT_NE(res, nullptr);

    EXPECT_EQ(res->name, std::string("PRO"));

    // Build ordered array: N, CA, CB, CG, CD, then remaining atoms
    // gemmi stores trimmed atom names (e.g. "N", "CA", "CB"), not padded
    std::vector<const gemmi::Atom*> ordered;
    const gemmi::Atom* n = nullptr;
    const gemmi::Atom* ca = nullptr;
    const gemmi::Atom* cb = nullptr;
    const gemmi::Atom* cg = nullptr;
    const gemmi::Atom* cd = nullptr;

    for (gemmi::Atom& atom : res->atoms) {
      std::string name(atom.name);
      if (name == "N") n = &atom;
      if (name == "CA") ca = &atom;
      if (name == "CB") cb = &atom;
      if (name == "CG") cg = &atom;
      if (name == "CD") cd = &atom;
    }

    EXPECT_TRUE(n && ca && cb && cg && cd);

    if (n && ca && cb && cg && cd) {
      ordered.push_back(n);
      ordered.push_back(ca);
      ordered.push_back(cb);
      ordered.push_back(cg);
      ordered.push_back(cd);
      for (gemmi::Atom& atom : res->atoms) {
        std::string name(atom.name);
        if (name == "N" || name == "CA" || name == "CB" || name == "CG" || name == "CD")
          continue;
        ordered.push_back(&atom);
      }
    }

    EXPECT_EQ(ordered.size(), 7u);

    std::vector<std::vector<int>> contacts =
      coot::util::get_contact_indices_for_PRO_residue_gemmi(ordered, &geom);

    EXPECT_EQ(contacts.size(), 7u);

    int total_contacts = 0;
    for (size_t i = 0; i < contacts.size(); i++) {
      total_contacts += contacts[i].size();
    }
    EXPECT_EQ(total_contacts, 7);

    // atom_0_contacts: 1,4
    ASSERT_EQ(contacts[0].size(), 2u);
    EXPECT_EQ(contacts[0][0], 1);
    EXPECT_EQ(contacts[0][1], 4);

    // atom_1_contacts: 5,2
    ASSERT_EQ(contacts[1].size(), 2u);
    EXPECT_EQ(contacts[1][0], 5);
    EXPECT_EQ(contacts[1][1], 2);

    // atom_2_contacts: 3
    ASSERT_EQ(contacts[2].size(), 1u);
    EXPECT_EQ(contacts[2][0], 3);

    // atom_3_contacts: 4
    ASSERT_EQ(contacts[3].size(), 1u);
    EXPECT_EQ(contacts[3][0], 4);

    // atom_4_contacts: (none)
    EXPECT_TRUE(contacts[4].empty());

    // atom_5_contacts: 6
    ASSERT_EQ(contacts[5].size(), 1u);
    EXPECT_EQ(contacts[5][0], 6);

    // atom_6_contacts: (none)
    EXPECT_TRUE(contacts[6].empty());
  }

  // Case 2: Non-PRO residue (ILE at //A/10)
  {
    std::string cid = "//A/10";
    gemmi::Residue* res = nullptr;
    for (gemmi::Chain& ch : model.chains) {
      if (ch.name == "A") {
        for (gemmi::Residue& r : ch.residues) {
          if (r.seqid.num.value == 10 && r.name == "ILE") {
            res = &r;
            break;
          }
        }
      }
    }
    ASSERT_NE(res, nullptr);
    EXPECT_EQ(res->name, std::string("ILE"));
  }

  // Case 3: Invalid CID
  {
    std::string cid = "//A/9999";
    gemmi::Residue* res = nullptr;
    for (gemmi::Chain& ch : model.chains) {
      if (ch.name == "A") {
        for (gemmi::Residue& r : ch.residues) {
          if (r.seqid.num.value == 9999) {
            res = &r;
            break;
          }
        }
      }
    }
    EXPECT_EQ(res, nullptr);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}