#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

TEST(OracleTestGemmi, make_header_metal_links_ng) {
  // Load PDB to verify chain/residue context
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  ASSERT_EQ(st.models.size(), 1u);

  gemmi::Model &model = st.models[0];
  ASSERT_GT(model.chains.size(), 0u);
  gemmi::Chain &chain = model.chains[0];

  // Find first ALA residue
  const gemmi::Residue *ala_res = nullptr;
  for (gemmi::Residue &r : chain.residues) {
    if (r.name == "ALA") {
      ala_res = &r;
      break;
    }
  }
  ASSERT_NE(ala_res, nullptr);

  std::string chain_id_str(chain.name);
  int seqnum = ala_res->seqid.num.value;
  std::string resname_str(ala_res->name);

  EXPECT_EQ(chain_id_str, "A");
  EXPECT_EQ(seqnum, 21);
  EXPECT_EQ(resname_str, "ALA");

  coot::protein_geometry pg;
  pg.init_standard();

  std::vector<link_info_t> links;

  // Case 1: ALA/O -> ZN/ZN (metal-O bond)
  {
    link_info_t link;
    link.chainID1  = chain_id_str;
    link.seqNum1   = seqnum;
    link.insCode1  = "";
    link.atName1   = "O";
    link.resName1  = resname_str;
    link.aloc1     = "";

    link.chainID2  = chain_id_str;
    link.seqNum2   = seqnum + 100;
    link.insCode2  = "";
    link.atName2   = "ZN";
    link.resName2  = "ZN";
    link.aloc2     = "";

    EXPECT_EQ(link.resName1, "ALA");
    EXPECT_EQ(link.atName1, "O");
    EXPECT_EQ(link.resName2, "ZN");
    EXPECT_EQ(link.atName2, "ZN");

    links.push_back(link);
    EXPECT_NO_THROW(coot::make_header_metal_links_ng_gemmi(links, pg));
  }

  // Count bonds after case 1 — original test expects 0
  {
    int bonds = coot::make_header_metal_links_ng_gemmi(links, pg);
    EXPECT_EQ(bonds, 0);
  }

  // Case 2: MET/SD -> ZN/ZN (metal-S bond)
  {
    link_info_t link;
    link.chainID1  = chain_id_str;
    link.seqNum1   = seqnum;
    link.insCode1  = "";
    link.atName1   = "SD";
    link.resName1  = "MET";
    link.aloc1     = "";

    link.chainID2  = chain_id_str;
    link.seqNum2   = seqnum + 200;
    link.insCode2  = "";
    link.atName2   = "ZN";
    link.resName2  = "ZN";
    link.aloc2     = "";

    EXPECT_EQ(link.resName1, "MET");
    EXPECT_EQ(link.atName1, "SD");
    EXPECT_EQ(link.resName2, "ZN");
    EXPECT_EQ(link.atName2, "ZN");

    links.push_back(link);
    EXPECT_NO_THROW(coot::make_header_metal_links_ng_gemmi(links, pg));
  }

  {
    int bonds = coot::make_header_metal_links_ng_gemmi(links, pg);
    EXPECT_EQ(bonds, 0);
  }

  // Case 3: ALA/O -> ALA/CB (non-metal, no bond added)
  {
    link_info_t link;
    link.chainID1  = chain_id_str;
    link.seqNum1   = seqnum;
    link.insCode1  = "";
    link.atName1   = "O";
    link.resName1  = resname_str;
    link.aloc1     = "";

    link.chainID2  = chain_id_str;
    link.seqNum2   = seqnum + 300;
    link.insCode2  = "";
    link.atName2   = "CB";
    link.resName2  = "ALA";
    link.aloc2     = "";

    EXPECT_EQ(link.resName1, "ALA");
    EXPECT_EQ(link.atName1, "O");
    EXPECT_EQ(link.resName2, "ALA");
    EXPECT_EQ(link.atName2, "CB");

    links.push_back(link);
    EXPECT_NO_THROW(coot::make_header_metal_links_ng_gemmi(links, pg));
  }

  {
    int bonds = coot::make_header_metal_links_ng_gemmi(links, pg);
    EXPECT_EQ(bonds, 0);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}