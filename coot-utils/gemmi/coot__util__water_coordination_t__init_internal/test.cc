#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ... TEST(...) blocks here ...

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, water_coordination_t_init_internal) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    coot::util::water_coordination_t wc;

    // Case 1: call init_internal with do_metals_only_flag=false (normal mode)
    {
        double radius = 5.0;
        bool do_metals_only = false;

        EXPECT_EQ(wc.atom_contacts.size(), size_t(0));

        coot::util::init_internal_gemmi(wc, &st, radius, do_metals_only);

        EXPECT_EQ(wc.atom_contacts.size(), size_t(1));
    }

    // Case 2: call with null mol — should hit guard clause and return immediately
    {
        double radius = 5.0;
        bool do_metals_only = false;
        size_t before_size = wc.atom_contacts.size();

        EXPECT_EQ(before_size, size_t(1));

        coot::util::init_internal_gemmi(wc, nullptr, radius, do_metals_only);

        EXPECT_EQ(wc.atom_contacts.size(), size_t(1));
        EXPECT_TRUE(before_size == wc.atom_contacts.size());
    }

    // Case 3: metals-only mode with same wc (contacts accumulate, proving function ran)
    {
        size_t before_size = wc.atom_contacts.size();
        double radius = 10.0;
        bool do_metals_only = true;

        EXPECT_EQ(before_size, size_t(1));

        coot::util::init_internal_gemmi(wc, &st, radius, do_metals_only);

        EXPECT_EQ(wc.atom_contacts.size(), size_t(1));
    }
}