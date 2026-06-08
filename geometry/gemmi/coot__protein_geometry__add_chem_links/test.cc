#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cifdoc.hpp>
#include "function.hh"

TEST(OracleTest, add_chem_links) {
    std::map<unsigned int, std::vector<coot::chem_link>> chem_link_map;

    // Case 1: Valid loop with 2 chem links
    {
        // Add 17 dummy entries with distinct chemical identities so each gets
        // its own hash bucket (hash is based on comp_id/mod_id, not link id)
        for (int i = 0; i < 17; i++) {
            std::string comp1 = "DUMM" + std::to_string(i);
            std::string comp2 = "DUMMYX" + std::to_string(i);
            std::string mod1 = "MODA" + std::to_string(i);
            std::string mod2 = "MODB" + std::to_string(i);
            coot::chem_link clink("dummy-id", comp1, mod1, comp1, comp2, mod2, comp2, "dummy-name");
            chem_link_map[clink.get_hash_code()].push_back(clink);
        }

        int before = chem_link_map.size();

        gemmi::cif::Loop loop;
        // Define the 8 columns by setting tags directly (public member)
        loop.tags = {"id", "comp_id_1", "mod_id_1", "group_comp_1",
                     "comp_id_2", "mod_id_2", "group_comp_2", "name"};

        // Row 0: first chem link
        loop.add_row({"link-1", "ALA", "N", "ALA", "GLY", "C", "GLY", "peptide"});

        // Row 1: second chem link
        loop.add_row({"link-2", "SER", "OG", "SER", "TYR", "OH", "TYR", "hydrogen-bond"});

        EXPECT_EQ(loop.length(), 2);
        EXPECT_EQ(before, 17);

        coot::add_chem_links_gemmi(loop, chem_link_map);

        int after = chem_link_map.size();
        EXPECT_EQ(after, 19);

        // Verify the two newly added chem links are present by scanning all entries
        bool found_link1 = false;
        bool found_link2 = false;
        for (auto const& kv : chem_link_map) {
            for (auto const& cl : kv.second) {
                if (cl.id == "link-1") {
                    found_link1 = true;
                    EXPECT_EQ(cl.chem_link_comp_id_1, "ALA");
                    EXPECT_EQ(cl.chem_link_mod_id_1, "N");
                    EXPECT_EQ(cl.chem_link_comp_id_2, "GLY");
                    EXPECT_EQ(cl.chem_link_mod_id_2, "C");
                    EXPECT_EQ(cl.chem_link_name, "peptide");
                }
                if (cl.id == "link-2") {
                    found_link2 = true;
                    EXPECT_EQ(cl.chem_link_comp_id_1, "SER");
                    EXPECT_EQ(cl.chem_link_mod_id_1, "OG");
                    EXPECT_EQ(cl.chem_link_comp_id_2, "TYR");
                    EXPECT_EQ(cl.chem_link_mod_id_2, "OH");
                    EXPECT_EQ(cl.chem_link_name, "hydrogen-bond");
                }
            }
        }
        EXPECT_TRUE(found_link1);
        EXPECT_TRUE(found_link2);
    }

    // Case 2: Empty loop — verifies no-op path
    {
        chem_link_map.clear();
        int before = chem_link_map.size();

        gemmi::cif::Loop loop;
        // Define columns but add no rows
        loop.tags = {"id", "comp_id_1", "mod_id_1", "group_comp_1",
                     "comp_id_2", "mod_id_2", "group_comp_2", "name"};

        EXPECT_EQ(loop.length(), 0);
        EXPECT_EQ(before, 0);

        coot::add_chem_links_gemmi(loop, chem_link_map);

        int after = chem_link_map.size();
        EXPECT_EQ(after, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
