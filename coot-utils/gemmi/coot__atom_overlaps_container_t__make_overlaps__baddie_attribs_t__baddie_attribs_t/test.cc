#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test the baddie_attribs_t constructor with gemmi types
TEST(BaddieAttribsTest, Constructor) {
    // Read a test structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    
    // Find two atoms to test
    gemmi::Atom* cr_at = nullptr;
    gemmi::Atom* n_at = nullptr;
    
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.name == "A") {
                    for (auto& atom : res.atoms) {
                        if (atom.name == " N  " && !cr_at) {
                            cr_at = &atom;
                        } else if (atom.name == " N  " && cr_at && cr_at != &atom) {
                            n_at = &atom;
                            break;
                        }
                    }
                }
                if (cr_at && n_at) break;
            }
            if (cr_at && n_at) break;
        }
        if (cr_at && n_at) break;
    }
    
    // Test the constructor
    if (cr_at && n_at) {
        coot::atom_overlaps_container_t::make_overlaps::baddie_attribs_t baddie(
            cr_at, n_at, 1.55f, 1.55f, 2.8f, 0.319354f, false, false
        );
        
        EXPECT_EQ(baddie.r_1, 1.55f);
        EXPECT_EQ(baddie.r_2, 1.55f);
        EXPECT_EQ(baddie.d, 2.8f);
        EXPECT_EQ(baddie.o, 0.319354f);
        EXPECT_EQ(baddie.is_hydrogen_bond, false);
        EXPECT_EQ(baddie.hydrogen_atom_is_first_atom, false);
        EXPECT_EQ(baddie.cr_at, cr_at);
        EXPECT_EQ(baddie.n_at, n_at);
    }
}

TEST(BaddieAttribsTest, Sorter) {
    // Test the static sorter function
    // sorter returns true if b1 should come before b2 (descending order by 'o')
    // b1.o = 0.5f, b2.o = 1.0f
    // Since we want descending order, b2 (with o=1.0) should come before b1 (with o=0.5)
    // So sorter(b1, b2) should return false (b1 should NOT come before b2)
    // And sorter(b2, b1) should return true (b2 should come before b1)
    coot::atom_overlaps_container_t::make_overlaps::baddie_attribs_t b1(nullptr, nullptr, 1.0f, 1.0f, 2.0f, 0.5f, false, false);
    coot::atom_overlaps_container_t::make_overlaps::baddie_attribs_t b2(nullptr, nullptr, 1.0f, 1.0f, 2.0f, 1.0f, false, false);
    
    // b2.o (1.0) > b1.o (0.5), so b2 should come first
    // sorter(b1, b2) = b2.o < b1.o = 1.0 < 0.5 = false
    EXPECT_FALSE(coot::atom_overlaps_container_t::make_overlaps::baddie_attribs_t::sorter(b1, b2));
    // sorter(b2, b1) = b1.o < b2.o = 0.5 < 1.0 = true
    EXPECT_TRUE(coot::atom_overlaps_container_t::make_overlaps::baddie_attribs_t::sorter(b2, b1));
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}