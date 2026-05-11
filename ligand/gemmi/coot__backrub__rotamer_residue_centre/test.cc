#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include "function.hh"

// Helper function to compare clipper::Coord_orth coordinates
static bool coord_equal(const clipper::Coord_orth &a, float x, float y, float z, float tol = 1e-4f) {
    return std::abs(a.x() - x) < tol && std::abs(a.y() - y) < tol && std::abs(a.z() - z) < tol;
}

TEST(OracleTest, BackrubRotamerResidueCentre) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(st.models.size(), 1) << "PDB loading failed";

    // Case 1: ALA 10 in chain A
    {
        gemmi::Model& model = st.models[0];
        gemmi::Chain* chain_a = nullptr;
        for (auto& c : model.chains) {
            if (c.name == "A") {
                chain_a = &c;
                break;
            }
        }
        ASSERT_NE(chain_a, nullptr) << "chain A not found";

        gemmi::Residue* res = nullptr;
        for (auto& r : chain_a->residues) {
            if (r.seqid.num.value == 10) {
                res = &r;
                break;
            }
        }
        ASSERT_NE(res, nullptr) << "get_residue failed for //A/10";

        gemmi::Residue* prev = nullptr;
        gemmi::Residue* next = nullptr;
        for (size_t i = 0; i < chain_a->residues.size(); ++i) {
            if (chain_a->residues[i].seqid.num.value == 9) prev = &chain_a->residues[i];
            if (chain_a->residues[i].seqid.num.value == 11) next = &chain_a->residues[i];
        }
        ASSERT_NE(prev, nullptr) << "neighbour prev missing for //A/10";
        ASSERT_NE(next, nullptr) << "neighbour next missing for //A/10";

        clipper::Coord_orth rc = coot::backrub::rotamer_residue_centre_gemmi(*res);

        // Expected values from oracle output
        EXPECT_NEAR(rc.x(), 27.2126, 1e-4);
        EXPECT_NEAR(rc.y(), 11.192, 1e-4);
        EXPECT_NEAR(rc.z(), 67.5274, 1e-4);
    }

    // Case 2: ASP 50 in chain A
    {
        gemmi::Model& model = st.models[0];
        gemmi::Chain* chain_a = nullptr;
        for (auto& c : model.chains) {
            if (c.name == "A") {
                chain_a = &c;
                break;
            }
        }
        ASSERT_NE(chain_a, nullptr) << "chain A not found";

        gemmi::Residue* res = nullptr;
        for (auto& r : chain_a->residues) {
            if (r.seqid.num.value == 50) {
                res = &r;
                break;
            }
        }
        ASSERT_NE(res, nullptr) << "get_residue failed for //A/50";

        gemmi::Residue* prev = nullptr;
        gemmi::Residue* next = nullptr;
        for (size_t i = 0; i < chain_a->residues.size(); ++i) {
            if (chain_a->residues[i].seqid.num.value == 49) prev = &chain_a->residues[i];
            if (chain_a->residues[i].seqid.num.value == 51) next = &chain_a->residues[i];
        }
        ASSERT_NE(prev, nullptr) << "neighbour prev missing for //A/50";
        ASSERT_NE(next, nullptr) << "neighbour next missing for //A/50";

        clipper::Coord_orth rc = coot::backrub::rotamer_residue_centre_gemmi(*res);

        // Expected values from oracle output
        EXPECT_NEAR(rc.x(), 6.15482, 1e-4);
        EXPECT_NEAR(rc.y(), -3.97091, 1e-4);
        EXPECT_NEAR(rc.z(), 56.4205, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}