#include <gtest/gtest.h>
#include <gemmi/cif.hpp>
#include <string>
#include <utility>
#include <vector>
#include "function.hh"

TEST(OracleTest, comp_chiral_gemmi) {
    // --- Case 1: SER — chiral residue with _chem_comp_chir loop ---
    {
        std::string cif_path = "/lmb/home/jdialpuri/Development/coot-dev/coot/monomers/s/SER.cif";
        int imol_enc = 0;

        auto cif_file = gemmi::cif::read_file(cif_path);

        // Find the _chem_comp_chir loop
        const gemmi::cif::Loop* chir_loop = nullptr;
        for (const auto& block : cif_file.blocks) {
            for (const auto& item : block.items) {
                if (item.type == gemmi::cif::ItemType::Loop) {
                    const auto& loop = item.loop;
                    if (!loop.tags.empty() && loop.tags[0] == "_chem_comp_chir.comp_id") {
                        chir_loop = &loop;
                        break;
                    }
                }
            }
            if (chir_loop) break;
        }

        EXPECT_NE(chir_loop, nullptr);
        if (chir_loop) {
            EXPECT_EQ(static_cast<int>(chir_loop->length()), 1);
        }

        auto result = coot::protein_geometry::comp_chiral_gemmi(*chir_loop, imol_enc);
        EXPECT_EQ(result.first, 1);
        EXPECT_EQ(result.second.size(), 1);
        EXPECT_EQ(result.second[0], "SER");
    }

    // --- Case 2: GLY — achiral residue, no _chem_comp_chir loop ---
    {
        std::string cif_path = "/lmb/home/jdialpuri/Development/coot-dev/coot/monomers/g/GLY.cif";
        int imol_enc = 0;

        auto cif_file = gemmi::cif::read_file(cif_path);

        // Find the _chem_comp_chir loop (should not exist)
        const gemmi::cif::Loop* chir_loop = nullptr;
        for (const auto& block : cif_file.blocks) {
            for (const auto& item : block.items) {
                if (item.type == gemmi::cif::ItemType::Loop) {
                    const auto& loop = item.loop;
                    if (!loop.tags.empty() && loop.tags[0] == "_chem_comp_chir.comp_id") {
                        chir_loop = &loop;
                        break;
                    }
                }
            }
            if (chir_loop) break;
        }

        EXPECT_EQ(chir_loop, nullptr);

        // When there is no chir loop, call comp_chiral_gemmi with an empty Loop
        // to replicate the original behavior where n_chiral == 0 and comp_ids_count == 0
        gemmi::cif::Loop empty_loop;
        auto result = coot::protein_geometry::comp_chiral_gemmi(empty_loop, imol_enc);
        EXPECT_EQ(result.first, 0);
        EXPECT_EQ(result.second.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}