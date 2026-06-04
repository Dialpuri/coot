#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, crankshaft) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];
    gemmi::Chain* chain = nullptr;
    for (auto& ch : model.chains) {
        if (ch.name == "A") { chain = &ch; break; }
    }
    EXPECT_TRUE(chain);

    // Find residues by sequence number
    gemmi::Residue* res_0 = nullptr;
    gemmi::Residue* res_1 = nullptr;
    gemmi::Residue* res_2 = nullptr;
    gemmi::Residue* res_3 = nullptr;

    for (auto& res : chain->residues) {
        if (res.seqid.num.value == 50 && res.name == "ARG") {
            res_0 = &res;
        }
        if (res.seqid.num.value == 51 && res.name == "GLU") {
            res_1 = &res;
        }
        if (res.seqid.num.value == 52 && res.name == "ILE") {
            res_2 = &res;
        }
        if (res.seqid.num.value == 53 && res.name == "SER") {
            res_3 = &res;
        }
    }

    EXPECT_STREQ(res_0->name.c_str(), "ARG");
    EXPECT_STREQ(res_1->name.c_str(), "GLU");
    EXPECT_STREQ(res_2->name.c_str(), "ILE");
    EXPECT_STREQ(res_3->name.c_str(), "SER");

    // Find CA atoms
    gemmi::Atom* ca0 = nullptr;
    gemmi::Atom* ca1 = nullptr;
    gemmi::Atom* ca2 = nullptr;
    gemmi::Atom* ca3 = nullptr;

    for (auto& at : res_0->atoms) {
        if (at.name == "CA") ca0 = &at;
    }
    for (auto& at : res_1->atoms) {
        if (at.name == "CA") ca1 = &at;
    }
    for (auto& at : res_2->atoms) {
        if (at.name == "CA") ca2 = &at;
    }
    for (auto& at : res_3->atoms) {
        if (at.name == "CA") ca3 = &at;
    }

    EXPECT_TRUE(ca0 && ca1 && ca2 && ca3);

    coot::crankshaft_set cs(*res_0, *res_1, *res_2, *res_3);

    EXPECT_EQ(cs.v.size(), 8u);
    EXPECT_NE(cs.ca_1, nullptr);
    EXPECT_NE(cs.ca_2, nullptr);

    float ang = M_PI / 6.0f;
    EXPECT_NEAR(ang, 0.523599f, 1e-4f);

    std::vector<float> orig_x(6, 0), orig_y(6, 0), orig_z(6, 0);
    int indices[] = {2, 3, 4, 5};
    for (int i = 0; i < 4; i++) {
        gemmi::Atom* at = cs.v[indices[i]];
        if (at) {
            orig_x[indices[i]] = at->pos.x;
            orig_y[indices[i]] = at->pos.y;
            orig_z[indices[i]] = at->pos.z;
        }
    }

    EXPECT_NEAR(orig_x[2], 12.914f, 1e-3f);
    EXPECT_NEAR(orig_y[2], -2.277f, 1e-3f);
    EXPECT_NEAR(orig_z[2], 56.102f, 1e-3f);

    cs.move_the_atoms_gemmi(ang);

    EXPECT_NEAR(cs.v[2]->pos.x, 12.9077f, 1e-3f);
    EXPECT_NEAR(cs.v[2]->pos.y, -2.00423f, 1e-3f);
    EXPECT_NEAR(cs.v[2]->pos.z, 56.0371f, 1e-3f);

    float dx = cs.v[2]->pos.x - orig_x[2];
    float dy = cs.v[2]->pos.y - orig_y[2];
    float dz = cs.v[2]->pos.z - orig_z[2];
    float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
    EXPECT_NEAR(dist, 0.280451f, 1e-4f);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}