#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, position_by_tetrahedron) {
    // Load PDB via gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];
    gemmi::Chain* chain = model.find_chain("A");

    // case: ALA residue 21 — has CA, C, N, CB
    {
        gemmi::Residue* ala21 = nullptr;
        for (gemmi::Residue& r : chain->residues) {
            if (r.seqid.num.value == 21 && r.name == "ALA") {
                ala21 = &r;
                break;
            }
        }

        EXPECT_NE(ala21, nullptr);

        const gemmi::Atom* at_ca = ala21->find_atom("CA", '*');
        const gemmi::Atom* at_c  = ala21->find_atom("C",  '*');
        const gemmi::Atom* at_n  = ala21->find_atom("N",  '*');
        const gemmi::Atom* at_cb = ala21->find_atom("CB", '*');

        EXPECT_NE(at_ca, nullptr);
        EXPECT_NE(at_c, nullptr);
        EXPECT_NE(at_n, nullptr);
        EXPECT_NE(at_cb, nullptr);

        clipper::Coord_orth pos = coot::reduce::position_by_tetrahedron_gemmi(at_ca, at_c, at_n, at_cb, 0.97);
        EXPECT_NEAR(pos.x(), 23.8997, 1e-4);
        EXPECT_NEAR(pos.y(), 6.60619, 1e-4);
        EXPECT_NEAR(pos.z(), 74.9858, 1e-4);
    }

    // case: VAL residue 7 — CB as central, CA/CG1/CG2 as neighbors
    {
        gemmi::Residue* val7 = nullptr;
        for (gemmi::Residue& r : chain->residues) {
            if (r.seqid.num.value == 7 && r.name == "VAL") {
                val7 = &r;
                break;
            }
        }

        EXPECT_NE(val7, nullptr);

        const gemmi::Atom* at_cb  = val7->find_atom("CB",  '*');
        const gemmi::Atom* at_ca  = val7->find_atom("CA",  '*');
        const gemmi::Atom* at_cg1 = val7->find_atom("CG1", '*');
        const gemmi::Atom* at_cg2 = val7->find_atom("CG2", '*');

        EXPECT_NE(at_cb, nullptr);
        EXPECT_NE(at_ca, nullptr);
        EXPECT_NE(at_cg1, nullptr);
        EXPECT_NE(at_cg2, nullptr);

        clipper::Coord_orth pos = coot::reduce::position_by_tetrahedron_gemmi(at_cb, at_ca, at_cg1, at_cg2, 1.09);
        EXPECT_NEAR(pos.x(), 26.7911, 1e-4);
        EXPECT_NEAR(pos.y(), 8.69056, 1e-4);
        EXPECT_NEAR(pos.z(), 75.1194, 1e-4);
    }

    // case: GLY residue 11 — no CB, so neighbor is null (edge case)
    {
        gemmi::Residue* gly11 = nullptr;
        for (gemmi::Residue& r : chain->residues) {
            if (r.seqid.num.value == 11 && r.name == "GLY") {
                gly11 = &r;
                break;
            }
        }

        EXPECT_NE(gly11, nullptr);

        const gemmi::Atom* at_ca = gly11->find_atom("CA", '*');
        const gemmi::Atom* at_c  = gly11->find_atom("C",  '*');
        const gemmi::Atom* at_n  = gly11->find_atom("N",  '*');
        const gemmi::Atom* at_cb = gly11->find_atom("CB", '*');

        EXPECT_NE(at_ca, nullptr);
        EXPECT_NE(at_c, nullptr);
        EXPECT_NE(at_n, nullptr);
        EXPECT_EQ(at_cb, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}