#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static std::string test_pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";

TEST(OracleTest, AveragePosition) {
    gemmi::Structure st = gemmi::read_pdb_file(test_pdb_path);
    ASSERT_FALSE(st.models.empty()) << "load failed";

    // case 1: valid residue with atoms
    {
        const gemmi::Model& model = st.models[0];
        bool found = false;
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& res : chain.residues) {
                if (chain.name == "A" && res.seqid.num.value == 10) {
                    found = true;
                    std::vector<clipper::Coord_orth> pts;
                    for (const gemmi::Atom& atom : res.atoms) {
                        clipper::Coord_orth pt(atom.pos.x, atom.pos.y, atom.pos.z);
                        pts.push_back(pt);
                    }

                    clipper::Coord_orth avg = coot::util::average_position_gemmi(pts);

                    EXPECT_NEAR(avg.x(), 27.2126, 1e-4);
                    EXPECT_NEAR(avg.y(), 11.192, 1e-4);
                    EXPECT_NEAR(avg.z(), 67.5274, 1e-4);
                    break;
                }
            }
            if (found) break;
        }
        ASSERT_TRUE(found) << "residue not found: //A/10";
    }

    // case 2: empty vector
    {
        std::vector<clipper::Coord_orth> pts;
        clipper::Coord_orth avg = coot::util::average_position_gemmi(pts);

        EXPECT_NEAR(avg.x(), 0, 1e-4);
        EXPECT_NEAR(avg.y(), 0, 1e-4);
        EXPECT_NEAR(avg.z(), 0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}