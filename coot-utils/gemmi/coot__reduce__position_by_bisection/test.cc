#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ReducePositionByBisectionGemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB file";

    const gemmi::Model& model = st.models[0];
    const gemmi::Chain* chain = nullptr;
    for (const auto& c : model.chains) {
        if (c.name == "A") {
            chain = &c;
            break;
        }
    }
    ASSERT_NE(chain, nullptr) << "Chain A not found";
    
    const gemmi::Residue* res = nullptr;
    for (const auto& r : chain->residues) {
        if (r.seqid.num.value == 3) {
            res = &r;
            break;
        }
    }
    ASSERT_NE(res, nullptr) << "Residue 3 not found";

    {
        const gemmi::Atom* at_1 = nullptr;
        const gemmi::Atom* at_2 = nullptr;
        const gemmi::Atom* at_3 = nullptr;
        
        for (const auto& atom : res->atoms) {
            if (atom.name == "N") at_1 = &atom;
            if (atom.name == "CA") at_2 = &atom;
            if (atom.name == "C") at_3 = &atom;
        }
        
        ASSERT_NE(at_1, nullptr);
        ASSERT_NE(at_2, nullptr);
        ASSERT_NE(at_3, nullptr);

        double bl = 0.96;
        clipper::Coord_orth result = coot::reduce::position_by_bisection_gemmi(at_1, at_2, at_3, bl);

        EXPECT_NEAR(result.x(), 15.7099, 1e-4);
        EXPECT_NEAR(result.y(), 0.782141, 1e-4);
        EXPECT_NEAR(result.z(), 78.7957, 1e-4);
    }

    {
        const gemmi::Atom* at_1_sc = nullptr;
        const gemmi::Atom* at_2_sc = nullptr;
        const gemmi::Atom* at_3_sc = nullptr;
        
        for (const auto& atom : res->atoms) {
            if (atom.name == "CG") at_1_sc = &atom;
            if (atom.name == "OD1") at_2_sc = &atom;
            if (atom.name == "ND2") at_3_sc = &atom;
        }
        
        ASSERT_NE(at_1_sc, nullptr);
        ASSERT_NE(at_2_sc, nullptr);
        ASSERT_NE(at_3_sc, nullptr);

        double bl_sc = 1.0;
        clipper::Coord_orth result_sc = coot::reduce::position_by_bisection_gemmi(at_1_sc, at_2_sc, at_3_sc, bl_sc);

        EXPECT_NEAR(result_sc.x(), 12.6032, 1e-4);
        EXPECT_NEAR(result_sc.y(), -0.111042, 1e-4);
        EXPECT_NEAR(result_sc.z(), 76.9498, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}