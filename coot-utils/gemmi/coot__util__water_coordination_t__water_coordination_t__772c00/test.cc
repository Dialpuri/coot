#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, water_coordination_t) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_GT(st.models.size(), 0u);

    // Case 1: Normal water coordination (all protein atoms), radius=5.0
    {
        double radius = 5.0;
        bool do_metals_only = false;

        coot::util::water_coordination_t wc = coot::util::water_coordination_t::water_coordination_t_gemmi(&st, radius, do_metals_only);

        auto contacts = wc.get_contacts();
        EXPECT_EQ(contacts.size(), 1u);

        auto metals = wc.metals();
        EXPECT_EQ(metals.size(), 0u);

        auto highly_coord = wc.get_highly_coordinated_waters(3, 5.0);
        EXPECT_EQ(highly_coord.size(), 1u);
    }

    // Case 2: Metals only, radius=5.0
    {
        double radius = 5.0;
        bool do_metals_only = true;

        coot::util::water_coordination_t wc = coot::util::water_coordination_t::water_coordination_t_gemmi(&st, radius, do_metals_only);

        auto contacts = wc.get_contacts();
        EXPECT_EQ(contacts.size(), 0u);

        auto metals = wc.metals();
        EXPECT_EQ(metals.size(), 0u);
    }

    // Case 3: Very large radius to catch more contacts
    {
        double radius = 10.0;
        bool do_metals_only = false;

        coot::util::water_coordination_t wc = coot::util::water_coordination_t::water_coordination_t_gemmi(&st, radius, do_metals_only);

        auto contacts = wc.get_contacts();
        EXPECT_EQ(contacts.size(), 1u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}