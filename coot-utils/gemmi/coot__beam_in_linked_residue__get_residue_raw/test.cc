#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include "function.hh"

// Helper to find a residue by chain ID and sequence number
static gemmi::Residue* find_residue_by_chain_seq(gemmi::Structure& st,
                                                  const std::string& chain_id,
                                                  int seqnum) {
  if (st.models.empty()) return nullptr;
  for (auto& chain : st.models[0].chains) {
    if (chain.name == chain_id) {
      for (auto& res : chain.residues) {
        if (res.seqid.num.value == seqnum) {
          return &res;
        }
      }
    }
  }
  return nullptr;
}

TEST(OracleTest, beam_in_linked_residue) {
    // Initialize protein geometry
    coot::protein_geometry geom;
    geom.init_standard();

    // Load carbohydrate residues (NAG and MAN)
    int nag_loaded = geom.try_dynamic_add("NAG", 1);
    int man_loaded = geom.try_dynamic_add("MAN", 1);
    EXPECT_EQ(nag_loaded, 1);
    EXPECT_EQ(man_loaded, 1);

    // Load PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    // Get a residue to use as residue_ref (chain A, seqnum 10)
    gemmi::Residue* residue_ref = find_residue_by_chain_seq(st, "A", 10);
    ASSERT_NE(residue_ref, nullptr);

    std::string ref_resname = residue_ref->name;
    int ref_seqnum = residue_ref->seqid.num.value;
    EXPECT_EQ(ref_resname, "ILE");
    EXPECT_EQ(ref_seqnum, 10);

    // Case 1: Valid sugar template (NAG - N-acetylglucosamine)
    {
        coot::beam_in_linked_residue beam(residue_ref, "BETA1-4", "NAG", &geom);

        // Check constructor-initialized public fields (non-trivial values)
        EXPECT_EQ(beam.link_type, "BETA1-4");
        EXPECT_EQ(beam.comp_id_new, "NAG");
        EXPECT_EQ(beam.comp_id_ref, ref_resname); // ILE (from residue_ref)
        EXPECT_NE(beam.residue_ref, nullptr);
        EXPECT_NE(beam.geom_p, nullptr);
        EXPECT_EQ(beam.template_res_ref, nullptr); // not set yet
        EXPECT_EQ(beam.template_res_mov, nullptr); // not set yet

        beam.setup_by_comp_id("NAG", "NAG");
        EXPECT_FALSE(beam.have_template);
        gemmi::Residue* result = beam.get_residue_raw_gemmi();
        EXPECT_EQ(result, nullptr);
    }

    // Case 2: Without setting up template - should return NULL
    {
        coot::beam_in_linked_residue beam(residue_ref, "ALPHA1-4", "NAG", &geom);
        EXPECT_EQ(beam.link_type, "ALPHA1-4");
        EXPECT_FALSE(beam.have_template);
        gemmi::Residue* result = beam.get_residue_raw_gemmi();
        EXPECT_EQ(result, nullptr);
    }

    // Case 3: Different sugar type (MAN - Mannose)
    {
        coot::beam_in_linked_residue beam(residue_ref, "BETA1-4", "MAN", &geom);
        EXPECT_EQ(beam.comp_id_new, "MAN");

        beam.setup_by_comp_id("MAN", "NAG");
        EXPECT_FALSE(beam.have_template);
        gemmi::Residue* result = beam.get_residue_raw_gemmi();
        EXPECT_EQ(result, nullptr);
    }

    // Case 4: Check setup_by_comp_id return value (bool) — should be false
    {
        coot::beam_in_linked_residue beam(residue_ref, "BETA1-4", "NAG", &geom);
        bool setup_result = beam.setup_by_comp_id("NAG", "NAG");
        EXPECT_FALSE(setup_result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}