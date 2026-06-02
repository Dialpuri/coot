#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include "function.hh"

TEST(OracleTest, Tree) {
    // Case 1: Tree with parent-child relationship (but SetAngles not called, so trivial output)
    {
        Tree tree;
        TreeVertex *v0 = new TreeVertex();
        TreeVertex *v1 = new TreeVertex();
        v0->id = 0;
        v1->id = 1;
        v0->coord = gemmi::Vec3(1.0, 2.0, 3.0);
        v1->coord = gemmi::Vec3(4.0, 5.0, 6.0);
        v1->parent = v0;
        tree.coords.push_back(v0);
        tree.coords.push_back(v1);

        EXPECT_EQ(tree.coords.size(), 2u);

        std::ostringstream oss;
        print_tree_gemmi(oss, tree);
        std::string output = oss.str();
        EXPECT_FALSE(output.empty());
        EXPECT_NE(output.find("Atom[0]"), std::string::npos);
        EXPECT_NE(output.find("Atom[1]"), std::string::npos);
        EXPECT_NE(output.find("1.000000"), std::string::npos);
        EXPECT_NE(output.find("has no parent"), std::string::npos);
        EXPECT_NE(output.find("Trees ...."), std::string::npos);
    }

    // Case 2: Tree with multiple independent vertices (no parents)
    {
        Tree tree;
        TreeVertex *a = new TreeVertex();
        TreeVertex *b = new TreeVertex();
        a->id = 0;
        b->id = 1;
        a->coord = gemmi::Vec3(0.0, 0.0, 0.0);
        b->coord = gemmi::Vec3(10.0, 20.0, 30.0);
        tree.coords.push_back(a);
        tree.coords.push_back(b);

        EXPECT_EQ(tree.coords.size(), 2u);

        std::ostringstream oss;
        print_tree_gemmi(oss, tree);
        std::string output = oss.str();
        EXPECT_FALSE(output.empty());
        EXPECT_NE(output.find("Atom[0]"), std::string::npos);
        EXPECT_NE(output.find("Atom[1]"), std::string::npos);
        EXPECT_NE(output.find("0.000000"), std::string::npos);
        EXPECT_NE(output.find("10.000000"), std::string::npos);
        EXPECT_NE(output.find("20.000000"), std::string::npos);
        EXPECT_NE(output.find("30.000000"), std::string::npos);
        EXPECT_NE(output.find("has no parent"), std::string::npos);
        EXPECT_NE(output.find("Trees ...."), std::string::npos);
    }

    // Case 3: Empty tree
    {
        Tree tree;
        EXPECT_EQ(tree.coords.size(), 0u);

        std::ostringstream oss;
        print_tree_gemmi(oss, tree);
        std::string output = oss.str();
        EXPECT_NE(output.find("Trees ...."), std::string::npos);
    }

    // Complementary case: Tree with 3-vertex chain to test non-trivial vertex count
    {
        Tree tree;
        TreeVertex *v0 = new TreeVertex();
        TreeVertex *v1 = new TreeVertex();
        TreeVertex *v2 = new TreeVertex();

        v0->id = 0;
        v1->id = 1;
        v2->id = 2;

        v0->coord = gemmi::Vec3(0.0, 0.0, 0.0);
        v1->coord = gemmi::Vec3(1.0, 0.0, 0.0);
        v2->coord = gemmi::Vec3(1.0, 1.0, 0.0);

        v1->parent = v0;
        v2->parent = v1;

        tree.coords.push_back(v0);
        tree.coords.push_back(v1);
        tree.coords.push_back(v2);

        EXPECT_EQ(tree.coords.size(), 3u);

        EXPECT_NE(v1->parent, nullptr);
        EXPECT_NE(v2->parent, nullptr);
        EXPECT_EQ(v0->parent, nullptr);

        std::ostringstream oss;
        print_tree_gemmi(oss, tree);
        std::string output = oss.str();
        EXPECT_FALSE(output.empty());
        EXPECT_NE(output.find("Atom[0]"), std::string::npos);
        EXPECT_NE(output.find("Atom[1]"), std::string::npos);
        EXPECT_NE(output.find("Atom[2]"), std::string::npos);
        EXPECT_NE(output.find("Trees ...."), std::string::npos);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}