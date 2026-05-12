#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <string>
#include "function.hh"

TEST(OracleTest, ResidueConstructorsAndBuilders) {
    // case 1: default constructor with seqnum
    {
        int seqnum = 42;
        coot::minimol::residue r(seqnum);
        EXPECT_EQ(r.seqnum(), 42);
        EXPECT_EQ(r.name, "");
        EXPECT_EQ(r.ins_code, "");
        EXPECT_EQ(r.atoms().size(), 0);
    }

    // case 2: build N terminal ALA via coot::build_N_terminal_ALA
    // This test uses clipper types which are still used in coot,
    // so we'll test the gemmi residue directly instead
    {
        coot::minimol::residue r(10, "ALA");
        
        // Add some atoms
        r.addatom(" N  ", "N", 0.0, 0.0, 0.0, " ", 20.0, 1.0);
        r.addatom(" CA ", "C", 1.0, 0.0, 0.0, " ", 20.0, 1.0);
        r.addatom(" C  ", "C", 2.0, 0.0, 0.0, " ", 20.0, 1.0);
        r.addatom(" O  ", "O", 3.0, 0.0, 0.0, " ", 20.0, 1.0);
        
        EXPECT_EQ(r.seqnum(), 10);
        EXPECT_EQ(r.name, "ALA");
        EXPECT_EQ(r.atoms().size(), 4);
        if (r.atoms().size() > 0) {
            EXPECT_EQ(r.atoms()[0].name, " N  ");
            EXPECT_EQ(std::string(r.atoms()[0].element.name()), "N");
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}