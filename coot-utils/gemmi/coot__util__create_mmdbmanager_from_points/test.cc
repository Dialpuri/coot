#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include "function.hh"

TEST(OracleTest, create_mmdbmanager_from_points) {
    // Case 1: typical — three points with a normal B-factor
    {
        std::vector<clipper::Coord_orth> pts;
        pts.push_back(clipper::Coord_orth(10.0f, 20.0f, 30.0f));
        pts.push_back(clipper::Coord_orth(40.0f, 50.0f, 60.0f));
        pts.push_back(clipper::Coord_orth(70.0f, 80.0f, 90.0f));
        float b_factor = 20.0f;

        gemmi::Structure st = coot::util::create_mmdbmanager_from_points_gemmi(pts, b_factor);

        EXPECT_EQ(st.models.size(), 1);
        EXPECT_EQ(st.models[0].chains.size(), 1);

        gemmi::Chain* chain = nullptr;
        for (auto& c : st.models[0].chains) {
            if (c.name == "A") { chain = &c; break; }
        }
        EXPECT_NE(chain, nullptr);
        if (chain) {
            EXPECT_EQ(chain->residues.size(), 3);
            EXPECT_EQ(gemmi::count_atom_sites(*chain), 3);
        }

        // MMDB GetResidue(1) = second residue (index 1, 0-based)
        gemmi::Residue& res = chain->residues[1];
        EXPECT_EQ(std::string(res.name), "ALA");
        EXPECT_EQ(res.seqid.num.value, 1);

        gemmi::Atom* atom = &res.atoms[0];
        if (atom) {
            EXPECT_EQ(atom->name, " CA ");
            std::string elem_str = std::string(" ") + atom->element.name();
            EXPECT_EQ(elem_str, " C");
        } else {
            FAIL() << "Expected atom to exist";
        }
    }

    // Case 2: edge — empty vector (no points)
    {
        std::vector<clipper::Coord_orth> pts;
        float b_factor = 0.0f;

        gemmi::Structure st = coot::util::create_mmdbmanager_from_points_gemmi(pts, b_factor);

        EXPECT_EQ(st.models.size(), 1);
        EXPECT_EQ(st.models[0].chains.size(), 1);

        gemmi::Chain* chain = nullptr;
        for (auto& c : st.models[0].chains) {
            if (c.name == "A") { chain = &c; break; }
        }
        EXPECT_NE(chain, nullptr);
        if (chain) {
            EXPECT_EQ(chain->residues.size(), 0);
        }
    }

    // Case 3: single point with high B-factor
    {
        std::vector<clipper::Coord_orth> pts;
        pts.push_back(clipper::Coord_orth(0.0f, 0.0f, 0.0f));
        float b_factor = 100.0f;

        gemmi::Structure st = coot::util::create_mmdbmanager_from_points_gemmi(pts, b_factor);

        EXPECT_EQ(st.models.size(), 1);

        gemmi::Chain* chain = nullptr;
        for (auto& c : st.models[0].chains) {
            if (c.name == "A") { chain = &c; break; }
        }
        EXPECT_NE(chain, nullptr);

        // MMDB GetResidue(0) = first residue (index 0)
        if (chain && !chain->residues.empty()) {
            gemmi::Residue& res = chain->residues[0];
            EXPECT_EQ(res.seqid.num.value, 0);
            EXPECT_EQ(std::string(res.name), "ALA");
        } else {
            FAIL() << "Expected residue to exist";
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}