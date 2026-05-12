#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(ArcInfoTest, Construction) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    const gemmi::Model& model = st.models[0];
    const gemmi::Chain* chain_a = nullptr;
    for (const auto& c : model.chains) {
        if (c.name == "A") {
            chain_a = &c;
            break;
        }
    }
    ASSERT_NE(chain_a, nullptr) << "Chain A not found";

    const gemmi::Residue* res_ile = nullptr;
    for (const auto& r : chain_a->residues) {
        if (r.seqid.num.value == 10) {
            res_ile = &r;
            break;
        }
    }
    ASSERT_NE(res_ile, nullptr) << "Residue ILE 10 not found";

    const gemmi::Atom* at_1 = nullptr;
    const gemmi::Atom* at_2 = nullptr;
    const gemmi::Atom* at_3 = nullptr;
    for (const auto& a : res_ile->atoms) {
        if (a.name == "CA") at_1 = &a;
        if (a.name == "CB") at_2 = &a;
        if (a.name == "CG1") at_3 = &a;
    }

    ASSERT_NE(at_1, nullptr) << "CA atom not found";
    ASSERT_NE(at_2, nullptr) << "CB atom not found";
    ASSERT_NE(at_3, nullptr) << "CG1 atom not found";

    {
        coot::arc_info_type arc(at_1, at_2, at_3);
        EXPECT_NEAR(arc.delta, 113.314, 1e-3);
        EXPECT_NEAR(arc.start_point.x, 27.686, 1e-3);
        EXPECT_NEAR(arc.start_point.y, 10.784, 1e-3);
        EXPECT_NEAR(arc.start_point.z, 67.801, 1e-3);
        EXPECT_NEAR(arc.normal.x, 0.239865, 1e-4);
        EXPECT_NEAR(arc.normal.y, 0.0235746, 1e-4);
        EXPECT_NEAR(arc.normal.z, -0.97052, 1e-4);
    }

    const gemmi::Residue* res_gly = nullptr;
    for (const auto& r : chain_a->residues) {
        if (r.seqid.num.value == 11) {
            res_gly = &r;
            break;
        }
    }
    ASSERT_NE(res_gly, nullptr) << "Residue GLY 11 not found";

    const gemmi::Atom* gly_at_1 = nullptr;
    const gemmi::Atom* gly_at_2 = nullptr;
    const gemmi::Atom* gly_at_3 = nullptr;
    for (const auto& a : res_gly->atoms) {
        if (a.name == "N") gly_at_1 = &a;
        if (a.name == "CA") gly_at_2 = &a;
        if (a.name == "C") gly_at_3 = &a;
    }

    ASSERT_NE(gly_at_1, nullptr) << "N atom not found";
    ASSERT_NE(gly_at_2, nullptr) << "CA atom not found";
    ASSERT_NE(gly_at_3, nullptr) << "C atom not found";

    {
        coot::arc_info_type arc(gly_at_1, gly_at_2, gly_at_3);
        EXPECT_NEAR(arc.delta, 112.627, 1e-3);
    }

    {
        EXPECT_ANY_THROW({
            coot::arc_info_type arc(at_1, at_2, nullptr);
        });
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}