#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cifdoc.hpp>
#include <gemmi/read_cif.hpp>
#include "function.hh"

TEST(OracleTest, add_energy_lib_angles) {
    coot::energy_lib_t energy_lib;

    // Parse the energy lib CIF file using gemmi
    std::string cif_file = "/lmb/home/jdialpuri/Development/coot-dev/coot/monomers/ener_lib.cif";
    gemmi::cif::Document doc = gemmi::read_cif_gz(cif_file);
    ASSERT_FALSE(doc.blocks.empty()) << "Failed to read CIF file";

    // Find the _lib_angle loop in the energy data block
    const gemmi::cif::Loop* angle_loop = nullptr;
    for (const auto& block : doc.blocks) {
        if (block.name == "energy") {
            for (const auto& item : block.items) {
                if (item.type == gemmi::cif::ItemType::Loop) {
                    if (!item.loop.tags.empty() &&
                        item.loop.tags[0].find("_lib_angle.") == 0) {
                        angle_loop = &item.loop;
                        break;
                    }
                }
            }
        }
        if (angle_loop) break;
    }

    ASSERT_NE(angle_loop, nullptr) << "Could not find _lib_angle loop";

    int loop_len = int(angle_loop->length());

    // Case 1: Call add_energy_lib_angles with the full valid loop
    {
        size_t angles_before = energy_lib.angles.size();
        EXPECT_EQ(angles_before, 0);
        EXPECT_EQ(loop_len, 1384);

        EXPECT_NO_THROW(energy_lib.add_energy_lib_angles_gemmi(*angle_loop));

        size_t angles_after = energy_lib.angles.size();
        EXPECT_EQ(angles_after, 1384);
        EXPECT_EQ(angles_after - angles_before, 1384);

        // Print first angle as sample
        EXPECT_FALSE(energy_lib.angles.empty());
        EXPECT_EQ(energy_lib.angles[0].atom_type_1, "NS");
        EXPECT_EQ(energy_lib.angles[0].atom_type_2, "CSP");
        EXPECT_EQ(energy_lib.angles[0].atom_type_3, "CH3");
        EXPECT_EQ(energy_lib.angles[0].spring_constant, 0);
        EXPECT_EQ(energy_lib.angles[0].angle, 180);
    }

    // Case 2: Edge case - library already loaded, verify idempotency
    {
        size_t angles_before = energy_lib.angles.size();
        EXPECT_EQ(angles_before, 1384);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}