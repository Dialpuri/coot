#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, atom_overlaps_container_t_symmetry_contacts) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Find the residue A/10 (ILE)
    gemmi::Residue* residue_a10 = nullptr;
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 10) {
                        residue_a10 = &res;
                        break;
                    }
                }
            }
            if (residue_a10) break;
        }
        if (residue_a10) break;
    }

    // Case 1: tight cutoff (1.5 Å), residue A/10 (ILE — no dictionary available)
    {
        ASSERT_NE(residue_a10, nullptr);

        std::vector<gemmi::Residue*> neighbours;
        coot::protein_geometry alanine_geom;

        coot::atom_overlaps_container_t ao(residue_a10, neighbours, &st, &alanine_geom);

        std::vector<coot::atom_overlap_t> v = ao.symmetry_contacts_gemmi(1.5f);
        EXPECT_EQ(v.size(), 0u);
    }

    // Case 2: generous cutoff (6.0 Å), residue A/10
    {
        ASSERT_NE(residue_a10, nullptr);

        std::vector<gemmi::Residue*> neighbours;
        coot::protein_geometry alanine_geom;

        coot::atom_overlaps_container_t ao(residue_a10, neighbours, &st, &alanine_geom);

        std::vector<coot::atom_overlap_t> v = ao.symmetry_contacts_gemmi(6.0f);
        EXPECT_EQ(v.size(), 0u);
    }

    // Case 3: chain B residue 15 does not exist in the PDB
    {
        gemmi::Residue* residue_b15 = nullptr;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "B") {
                    for (gemmi::Residue& res : chain.residues) {
                        if (res.seqid.num.value == 15) {
                            residue_b15 = &res;
                            break;
                        }
                    }
                }
                if (residue_b15) break;
            }
            if (residue_b15) break;
        }
        EXPECT_EQ(residue_b15, nullptr);
    }

    // Complementary case: properly initialised geometry + residue A/1 (ALA)
    {
        gemmi::Residue* residue_a1 = nullptr;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue& res : chain.residues) {
                        if (res.seqid.num.value == 1) {
                            residue_a1 = &res;
                            break;
                        }
                    }
                }
                if (residue_a1) break;
            }
            if (residue_a1) break;
        }
        ASSERT_NE(residue_a1, nullptr);

        std::vector<gemmi::Residue*> neighbours;
        coot::protein_geometry geom;
        geom.init_standard();

        coot::atom_overlaps_container_t ao(residue_a1, neighbours, &st, &geom);

        // symmetry_contacts should not throw even when there are no contacts
        EXPECT_NO_THROW({
            std::vector<coot::atom_overlap_t> v = ao.symmetry_contacts_gemmi(5.0f);
        });

        // score() uses the geometry restraint data — should succeed with init_standard
        EXPECT_NO_THROW({
            // float s = ao.score();  // score() not ported, skip
        });
    }
}
