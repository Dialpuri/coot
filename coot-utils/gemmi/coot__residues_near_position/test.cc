#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ResiduesNearPosition) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: coordinate within the protein — should find nearby residues
    clipper::Coord_orth pt1(12.0, 22.0, 31.0);
    double radius1 = 5.0;

    std::cout << "INPUT  coord: " << pt1.x() << " " << pt1.y() << " " << pt1.z() << std::endl;
    std::cout << "INPUT  radius: " << radius1 << std::endl;

    std::vector<gemmi::Residue*> v1 = coot::residues_near_position(pt1, st, radius1);
    std::cout << "OUTPUT count: " << v1.size() << std::endl;
    EXPECT_EQ(v1.size(), 1u);

    if (v1.size() > 0) {
        gemmi::Residue* res = v1[0];
        int seqnum = res->seqid.num.value;
        std::string chain_id;
        std::string resname = res->name;
        
        // Search for chain containing this residue
        for (auto& chain : st.models[0].chains) {
            for (auto& r : chain.residues) {
                if (&r == res) {
                    chain_id = chain.name;
                    break;
                }
            }
            if (!chain_id.empty()) break;
        }
        
        std::cout << "OUTPUT residue[0]: " << chain_id << " " << seqnum << " " << resname << std::endl;
        EXPECT_EQ(chain_id, "A");
        EXPECT_EQ(seqnum, 237);
        EXPECT_EQ(resname, "LYS");
    }

    // Case 2: coordinate far outside the protein — should return empty vector
    {
        clipper::Coord_orth pt2(1000.0, 1000.0, 1000.0);
        double radius2 = 0.5;
        std::cout << "INPUT  coord: " << pt2.x() << " " << pt2.y() << " " << pt2.z() << std::endl;
        std::cout << "INPUT  radius: " << radius2 << std::endl;
        std::vector<gemmi::Residue*> v2 = coot::residues_near_position(pt2, st, radius2);
        std::cout << "OUTPUT count: " << v2.size() << std::endl;
        EXPECT_EQ(v2.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}