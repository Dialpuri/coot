#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "function.hh"

TEST(OracleTest, PositionPairByBisection) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    gemmi::Residue* res = nullptr;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name != "A") continue;
            for (auto& r : chain.residues) {
                if (r.name == "GLU" && r.seqid.num.value == 2) {
                    res = &r;
                    break;
                }
            }
            if (res) break;
        }
        if (res) break;
    }
    ASSERT_NE(res, nullptr) << "Residue not found";
    
    const gemmi::Atom* at_1 = nullptr;
    const gemmi::Atom* at_2 = nullptr;
    const gemmi::Atom* at_3 = nullptr;
    
    for (const auto& atom : res->atoms) {
        if (atom.name == "CA") at_1 = &atom;
        else if (atom.name == "CB") at_2 = &atom;
        else if (atom.name == "CG") at_3 = &atom;
    }
    
    ASSERT_NE(at_1, nullptr) << "CA atom missing";
    ASSERT_NE(at_2, nullptr) << "CB atom missing";
    ASSERT_NE(at_3, nullptr) << "CG atom missing";
    
    double bl = 1.0;
    double alpha = 109.5 * M_PI / 180.0;
    
    std::pair<clipper::Coord_orth, clipper::Coord_orth> result = 
        coot::reduce::position_pair_by_bisection_gemmi(at_1, at_2, at_3, bl, alpha);
    
    EXPECT_NEAR(result.first.x(), 12.5475, 1e-4);
    EXPECT_NEAR(result.first.y(), 5.2838, 1e-4);
    EXPECT_NEAR(result.first.z(), 80.1273, 1e-4);
    EXPECT_NEAR(result.second.x(), 11.5465, 1e-4);
    EXPECT_NEAR(result.second.y(), 4.9369, 1e-4);
    EXPECT_NEAR(result.second.z(), 78.8843, 1e-4);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}