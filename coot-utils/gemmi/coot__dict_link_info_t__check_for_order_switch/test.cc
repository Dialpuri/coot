#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <string>
#include "geometry/protein-geometry.hh"

// Helper: find a residue in a specific chain by sequence number
static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seq_num) {
  for (auto& chain : model.chains) {
    if (chain.name == chain_id) {
      for (auto& res : chain.residues) {
        if (res.seqid.num.value == seq_num) {
          return &res;
        }
      }
    }
  }
  return nullptr;
}

TEST(OracleTest, check_for_order_switch) {
    // Load the PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_GE(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];

    // Init protein geometry with standard dictionary
    coot::protein_geometry geom;
    geom.init_standard();

    // Get residues for testing - oracle asked for 10/11 but got 11(GLY)/12(GLU)
    gemmi::Residue* res10 = find_residue(model, "A", 11);
    gemmi::Residue* res11 = find_residue(model, "A", 12);

    ASSERT_NE(res10, nullptr);
    ASSERT_NE(res11, nullptr);

    // Verify oracle's actual residues (oracle asked for 10/11 but got 11/GLY and 12/GLU)
    EXPECT_EQ(res10->seqid.num.value, 11);
    EXPECT_EQ(res11->seqid.num.value, 12);

    // Construct dict_link_info_t
    coot::dict_link_info_t info;

    // Case 1: Standard peptide bond TRANS between GLY11 and GLU12 (normal order)
    {
        std::string link_type = "TRANS";
        bool result = info.check_for_order_switch_gemmi(res10, res11, link_type, geom);
        EXPECT_FALSE(result);
    }

    // Case 2: Same residues but reversed order
    {
        std::string link_type = "TRANS";
        bool result = info.check_for_order_switch_gemmi(res11, res10, link_type, geom);
        EXPECT_FALSE(result);
    }

    // Case 3: Non-existent link type
    {
        std::string link_type = "NONEXISTENT_LINK";
        bool result = info.check_for_order_switch_gemmi(res10, res11, link_type, geom);
        EXPECT_FALSE(result);
    }

    // Case 4 (complementary): Try with residues that might produce a true order switch.
    // We need a case where matching_chem_links returns links but none match the
    // specific order of ref/new. Let's try with HOH (water) residues that may
    // have directional links.
    {
        // Try HOH residues - water molecules often have directional hydrogen bonds
        gemmi::Residue* hoh1 = find_residue(model, "A", 113);
        gemmi::Residue* hoh2 = find_residue(model, "A", 114);
        if (hoh1 && hoh2) {
            std::string link_type = "HOH";
            bool result = info.check_for_order_switch_gemmi(hoh1, hoh2, link_type, geom);
            // Just observe the result - may be true or false depending on dict
            (void)result;
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}