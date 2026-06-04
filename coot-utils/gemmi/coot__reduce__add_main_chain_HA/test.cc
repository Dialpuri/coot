#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static gemmi::Residue * find_residue_by_cid(gemmi::Model &model, const std::string &cid, std::string &out_chain_id) {
  // CID format: "//A/10" -> chain "A", seqnum 10
  // Parse chain and seqnum from CID
  std::string chain_id;
  int seq_num = -1;
  
  size_t first_slash = cid.find('/');
  size_t second_slash = cid.find('/', first_slash + 1);
  size_t third_slash = cid.find('/', second_slash + 1);
  
  if (second_slash != std::string::npos && third_slash != std::string::npos) {
    chain_id = cid.substr(second_slash + 1, third_slash - second_slash - 1);
    seq_num = std::stoi(cid.substr(third_slash + 1));
  }
  
  for (auto &chain : model.chains) {
    if (chain.name == chain_id) {
      for (auto &res : chain.residues) {
        if (res.seqid.num.value == seq_num) {
          out_chain_id = chain.name;
          return &res;
        }
      }
    }
  }
  return nullptr;
}

static std::string residue_get_res_name(const gemmi::Residue &res) {
  return res.name;
}

static int residue_get_seq_num(const gemmi::Residue &res) {
  return res.seqid.num.value;
}

static int residue_get_number_of_atoms(const gemmi::Residue &res) {
  return static_cast<int>(res.atoms.size());
}

static const gemmi::Atom * residue_get_atom(const gemmi::Residue &res, const std::string &atom_name) {
  for (const auto &atom : res.atoms) {
    if (atom.name == atom_name) {
      return &atom;
    }
  }
  return nullptr;
}

TEST(OracleTest, add_main_chain_HA) {
    // Setup: load PDB
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
        gemmi::Model &model = st.models[0];

        // Case 1: Non-glycine residue with HA (ILE at //A/10)
        {
            std::string cid = "//A/10";
            std::string chain_id;
            gemmi::Residue *res = find_residue_by_cid(model, cid, chain_id);
            EXPECT_NE(res, nullptr);
            EXPECT_STREQ(residue_get_res_name(*res).c_str(), "ILE");
            EXPECT_STREQ(chain_id.c_str(), "A");
            EXPECT_EQ(residue_get_seq_num(*res), 10);

            int atoms_before = residue_get_number_of_atoms(*res);
            EXPECT_EQ(atoms_before, 19);

            coot::reduce::add_main_chain_HA_gemmi(*res);

            int atoms_after = residue_get_number_of_atoms(*res);
            EXPECT_EQ(atoms_after, 19);
            EXPECT_EQ((atoms_after - atoms_before), 0);

            // Check if HA was added/found
            const gemmi::Atom *ha = residue_get_atom(*res, " HA ");
            EXPECT_NE(ha, nullptr);
            EXPECT_NEAR(ha->pos.x, 27.9124, 1e-4);
            EXPECT_NEAR(ha->pos.y, 12.7612, 1e-4);
            EXPECT_NEAR(ha->pos.z, 68.0598, 1e-4);
        }

        // Case 2: GLU residue (also non-glycine) at //A/2
        {
            std::string cid2 = "//A/2";
            std::string chain_id2;
            gemmi::Residue *res2 = find_residue_by_cid(model, cid2, chain_id2);
            EXPECT_NE(res2, nullptr);
            EXPECT_STREQ(residue_get_res_name(*res2).c_str(), "GLU");

            int atoms_before2 = residue_get_number_of_atoms(*res2);
            EXPECT_EQ(atoms_before2, 15);

            coot::reduce::add_main_chain_HA_gemmi(*res2);

            int atoms_after2 = residue_get_number_of_atoms(*res2);
            EXPECT_EQ(atoms_after2, 15);
            EXPECT_EQ((atoms_after2 - atoms_before2), 0);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}