#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

struct test_links_container {
    std::vector<gemmi::Connection> connections;

    void fill_links_gemmi(const gemmi::Structure* st) {
        connections.clear();
        if (st) {
            for (const auto& conn : st->connections) {
                connections.push_back(conn);
            }
        }
    }
};

TEST(OracleTest, fill_links_gemmi) {
    // Case: setup - load PDB
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    test_links_container container;

    // Case: initial_links_size
    EXPECT_EQ(container.connections.size(), 0);

    // Case: call fill_links_gemmi()
    container.fill_links_gemmi(&st);

    // Case: final_links_size
    EXPECT_EQ(container.connections.size(), 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}