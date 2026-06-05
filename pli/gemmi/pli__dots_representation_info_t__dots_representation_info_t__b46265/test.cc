#include <gtest/gtest.h>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, dots_representation_info_t) {
    gemmi::cif::Document doc = gemmi::cif::read_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    gemmi::Structure st = gemmi::make_structure(std::move(doc));

    // Case 1: default constructor (no mol) — checks empty init
    {
        pli::dots_representation_info_t dots_default;
        EXPECT_TRUE(dots_default.is_open_p());
        EXPECT_EQ(dots_default.get_name(), "");
    }

    // Case 2: constructor with gemmi::Structure* — the target overload
    {
        pli::dots_representation_info_t dots = pli::dots_representation_info_t::dots_representation_info_t_gemmi(&st);
        EXPECT_TRUE(dots.is_open_p());
        EXPECT_EQ(dots.get_name(), "");
    }

    // Case 3: constructor with string name
    {
        std::string test_name = "MyDots";
        pli::dots_representation_info_t dots_named(test_name);
        EXPECT_TRUE(dots_named.is_open_p());
        EXPECT_EQ(dots_named.get_name(), "MyDots");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}