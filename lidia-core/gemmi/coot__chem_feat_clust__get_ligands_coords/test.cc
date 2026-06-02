#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_ligands_coords) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // Find first residue in chain A
    const gemmi::Residue* first_res = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            if (!chain.residues.empty()) {
                first_res = &chain.residues[0];
            }
            break;
        }
    }

    ASSERT_NE(first_res, nullptr);

    // Case 1: Empty ligands vector — should return empty coord list
    {
        coot::chem_feat_clust::chem_feat_clust_t clust;
        std::vector<clipper::Coord_orth> coords = coot::chem_feat_clust::get_ligands_coords_gemmi(clust);
        EXPECT_EQ(clust.ligands.size(), 0u);   // INPUT ligands_count: 0
        EXPECT_EQ(coords.size(), 0u);          // OUTPUT coords_count: 0
    }

    // Case 2: One ligand with null residue — should return empty coord list
    {
        coot::chem_feat_clust::chem_feat_clust_t clust;

        coot::chem_feat_clust::ligand_t lig;
        lig.residue = nullptr;

        clust.ligands.push_back(lig);

        std::vector<clipper::Coord_orth> coords = coot::chem_feat_clust::get_ligands_coords_gemmi(clust);
        EXPECT_EQ(clust.ligands.size(), 1u);           // INPUT ligands_count: 1
        EXPECT_EQ(lig.residue, nullptr);               // INPUT residue_ptr: nullptr
        EXPECT_EQ(coords.size(), 0u);                  // OUTPUT coords_count: 0
    }

    // Case 3: One ligand with a valid residue — should return atom coords
    {
        coot::chem_feat_clust::chem_feat_clust_t clust;

        coot::chem_feat_clust::ligand_t lig;
        lig.residue = const_cast<gemmi::Residue*>(first_res);

        clust.ligands.push_back(lig);

        std::vector<clipper::Coord_orth> coords = coot::chem_feat_clust::get_ligands_coords_gemmi(clust);

        EXPECT_EQ(clust.ligands.size(), 1u);                                  // INPUT ligands_count: 1
        EXPECT_EQ(first_res->seqid.num.value, 1);                             // INPUT residue_seqnum: 1
        EXPECT_EQ(first_res->name, "MET");                                    // INPUT residue_name: MET
        EXPECT_TRUE(clust.ligands[0].residue);                                // INPUT residue_found: true
        EXPECT_EQ(coords.size(), 8u);                                         // OUTPUT coords_count: 8
        EXPECT_NEAR(coords[0][0], 12.334, 1e-4);                              // OUTPUT first_coord_x: 12.334
        EXPECT_NEAR(coords[0][1], 2.772, 1e-4);                               // OUTPUT first_coord_y: 2.772
        EXPECT_NEAR(coords[0][2], 74.951, 1e-4);                              // OUTPUT first_coord_z: 74.951
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}