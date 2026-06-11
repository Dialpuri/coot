#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include "function.hh"

#include "coot-utils/coot-map-utils.hh"

#include <cstdlib>

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int r = RUN_ALL_TESTS();
    std::fflush(nullptr);
    _exit(r);
}

TEST(OracleTest, make_gaussian_atom_map_for_mask) {
    // Load MTZ for reference map (provides cell/gridding)
    // Use molecules_container_t to get correct sampling rate
    clipper::Xmap<float> xmap;
    bool ok = coot::util::map_fill_from_mtz(
        &xmap,
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
        "FWT", "PHWT", "", false, 1.8f);
    ASSERT_TRUE(ok);
    const clipper::Xmap<float> &xmap_ref = xmap;

    // Load CIF for molecule using gemmi
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    gemmi::Model &model = st.models[0];

    float sigma = 1.5f;
    float box_radius = 5.0f;

    // Case 1: Known residue ILE 10 (8 atoms) — using CID
    {
        std::string cid = "//A/10";

        // Count atoms in the selection for verification
        int n_atoms = 0;
        for (gemmi::Chain &chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue &residue : chain.residues) {
                    if (residue.seqid.num.value == 10) {
                        n_atoms = static_cast<int>(residue.atoms.size());
                    }
                }
            }
        }
        EXPECT_EQ(n_atoms, 8);

        clipper::Xmap<float> result = coot::util::make_gaussian_atom_map_for_mask_gemmi(
            xmap_ref, model, cid, sigma, box_radius);

        float map_sum = 0;
        int count = 0;
        for (clipper::Xmap_base::Map_reference_coord iu = result.first_coord(); !iu.last(); iu.next()) {
            map_sum += result[iu];
            count++;
        }
        EXPECT_NEAR(map_sum, 1223.65f, 0.01f);
        EXPECT_EQ(count, 559872);
    }

    // Case 2: LZA ligand — try with full atom CID
    {
        std::string cid = "//A/1299";

        // Count atoms for verification
        int n_atoms = 0;
        for (gemmi::Chain &chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue &residue : chain.residues) {
                    if (residue.seqid.num.value == 1299) {
                        n_atoms = static_cast<int>(residue.atoms.size());
                    }
                }
            }
        }
        EXPECT_EQ(n_atoms, 25);

        clipper::Xmap<float> result = coot::util::make_gaussian_atom_map_for_mask_gemmi(
            xmap_ref, model, cid, sigma, box_radius);

        float map_sum = 0;
        for (clipper::Xmap_base::Map_reference_coord iu = result.first_coord(); !iu.last(); iu.next()) {
            map_sum += result[iu];
        }
        EXPECT_NEAR(map_sum, 3823.9f, 0.01f);
    }

    // Case 3: Non-existent residue
    {
        std::string cid = "//A/99999";

        // Count atoms for verification (should be 0)
        int n_atoms = 0;
        for (gemmi::Chain &chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue &residue : chain.residues) {
                    if (residue.seqid.num.value == 99999) {
                        n_atoms = static_cast<int>(residue.atoms.size());
                    }
                }
            }
        }
        EXPECT_EQ(n_atoms, 0);

        clipper::Xmap<float> result = coot::util::make_gaussian_atom_map_for_mask_gemmi(
            xmap_ref, model, cid, sigma, box_radius);

        float map_sum = 0;
        for (clipper::Xmap_base::Map_reference_coord iu = result.first_coord(); !iu.last(); iu.next()) {
            map_sum += result[iu];
        }
        EXPECT_NEAR(map_sum, 0.0f, 1e-4);
    }
}
