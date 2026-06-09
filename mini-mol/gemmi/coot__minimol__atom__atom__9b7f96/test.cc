#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static gemmi::Atom* find_atom_by_cid(gemmi::Model& model, std::string cid) {
  // Parse CID: "//A/1/CA" -> chain="A", resno=1, atom_name="CA"
  // Format is //chain/resno/atomname
  std::string chain_id;
  int resno = 0;
  std::string atom_name;

  size_t pos = 2; // skip "//"
  size_t slash = cid.find('/', pos);
  if (slash == std::string::npos) return nullptr;
  chain_id = cid.substr(pos, slash - pos);
  pos = slash + 1;

  slash = cid.find('/', pos);
  if (slash == std::string::npos) return nullptr;
  resno = std::stoi(cid.substr(pos, slash - pos));
  pos = slash + 1;
  atom_name = cid.substr(pos);

  for (gemmi::Chain& chain : model.chains) {
    if (chain.name != chain_id) continue;
    for (gemmi::Residue& res : chain.residues) {
      if (res.seqid.num.value != resno) continue;
      for (gemmi::Atom& at : res.atoms) {
        if (at.name == atom_name) return &at;
      }
    }
  }
  return nullptr;
}

TEST(OracleTest, cid_to_atom) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  gemmi::Model& model = st.models[0];

  // Case 1: valid atom — //A/1/CA
  {
    std::string cid = "//A/1/CA";
    gemmi::Atom* at = find_atom_by_cid(model, cid);
    ASSERT_NE(at, nullptr);
    coot::minimol::atom atom = coot::minimol::atom_gemmi(*at);
    EXPECT_EQ(atom.name, " CA ");
    EXPECT_EQ(atom.element, " C");
    EXPECT_NEAR(atom.pos.x(), 13.401, 1e-4);
    EXPECT_NEAR(atom.pos.y(), 3.78, 1e-4);
    EXPECT_NEAR(atom.pos.z(), 74.703, 1e-4);
    EXPECT_EQ(atom.altLoc, "");
    EXPECT_NEAR(atom.occupancy, 1.0, 1e-4);
    EXPECT_NEAR(atom.temperature_factor, 69.04, 1e-4);
    EXPECT_EQ(atom.int_user_data, -1);
  }

  // Case 2: valid atom — //A/50/CA
  {
    std::string cid = "//A/50/CA";
    gemmi::Atom* at = find_atom_by_cid(model, cid);
    ASSERT_NE(at, nullptr);
    coot::minimol::atom atom = coot::minimol::atom_gemmi(*at);
    EXPECT_EQ(atom.name, " CA ");
    EXPECT_EQ(atom.element, " C");
    EXPECT_NEAR(atom.pos.x(), 8.709, 1e-4);
    EXPECT_NEAR(atom.pos.y(), -3.609, 1e-4);
    EXPECT_NEAR(atom.pos.z(), 56.949, 1e-4);
    EXPECT_EQ(atom.altLoc, "");
    EXPECT_NEAR(atom.occupancy, 1.0, 1e-4);
    EXPECT_NEAR(atom.temperature_factor, 58.72, 1e-4);
    EXPECT_EQ(atom.int_user_data, -1);
  }

  // Case 3: invalid CID — non-existent atom
  {
    std::string cid = "//A/9999/CA";
    gemmi::Atom* at = find_atom_by_cid(model, cid);
    EXPECT_EQ(at, nullptr);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
