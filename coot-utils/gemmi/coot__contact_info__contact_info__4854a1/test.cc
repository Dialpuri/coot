#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ... TEST(...) blocks here ...

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, contact_info) {
    // Case 1: Valid contacts with 3 entries
    {
        std::vector<std::pair<int, int>> contacts_vec;
        contacts_vec.push_back({10, 15});
        contacts_vec.push_back({20, 25});
        contacts_vec.push_back({30, 35});

        coot::contact_info ci = coot::contact_info::contact_info_gemmi(contacts_vec);
        EXPECT_EQ(ci.n_contacts(), 3);
        EXPECT_EQ(ci.contacts[0].id1, 10);
        EXPECT_EQ(ci.contacts[0].id2, 15);
        EXPECT_EQ(ci.contacts[1].id1, 20);
        EXPECT_EQ(ci.contacts[1].id2, 25);
        EXPECT_EQ(ci.contacts[2].id1, 30);
        EXPECT_EQ(ci.contacts[2].id2, 35);
    }

    // Case 2: Zero contacts (edge case)
    {
        std::vector<std::pair<int, int>> contacts_vec;

        coot::contact_info ci = coot::contact_info::contact_info_gemmi(contacts_vec);
        EXPECT_EQ(ci.n_contacts(), 0);
    }

    // Case 3: Single contact
    {
        std::vector<std::pair<int, int>> contacts_vec;
        contacts_vec.push_back({100, 200});

        coot::contact_info ci = coot::contact_info::contact_info_gemmi(contacts_vec);
        EXPECT_EQ(ci.n_contacts(), 1);
        EXPECT_EQ(ci.contacts[0].id1, 100);
        EXPECT_EQ(ci.contacts[0].id2, 200);
    }
}