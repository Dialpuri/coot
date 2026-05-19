#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "function.hh"

// Helper to extract atoms from gemmi Residue as clipper::Coord_orth
static std::vector<clipper::Coord_orth>
get_atom_coords(const gemmi::Residue& res) {
    std::vector<clipper::Coord_orth> coords;
    coords.reserve(res.atoms.size());
    for (const auto& atom : res.atoms) {
        coords.emplace_back(atom.pos.x, atom.pos.y, atom.pos.z);
    }
    return coords;
}

// Helper to compute centroid of atoms
static clipper::Coord_orth
compute_centroid(const std::vector<clipper::Coord_orth>& coords) {
    float sum_x = 0.0f, sum_y = 0.0f, sum_z = 0.0f;
    if (!coords.empty()) {
        for (const auto& c : coords) {
            sum_x += c.x();
            sum_y += c.y();
            sum_z += c.z();
        }
        float inv_n = 1.0f / static_cast<float>(coords.size());
        return clipper::Coord_orth(sum_x * inv_n, sum_y * inv_n, sum_z * inv_n);
    }
    return clipper::Coord_orth(0, 0, 0);
}

TEST(OracleTest, backrub_residue_radius) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];
    bool found_chain = false;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            found_chain = true;

            // Case 1: residue 10 (GLY)
            {
                gemmi::Residue* this_res = nullptr;
                for (gemmi::Residue& r : chain.residues) {
                    if (r.seqid.num.value == 10) {
                        this_res = &r;
                        break;
                    }
                }
                ASSERT_NE(this_res, nullptr);

                std::vector<clipper::Coord_orth> residue_atoms = get_atom_coords(*this_res);
                clipper::Coord_orth rc = compute_centroid(residue_atoms);

                // Call our gemmi port
                float r = coot::backrub::residue_radius_gemmi(residue_atoms, rc);
                EXPECT_NEAR(r, 2.89078f, 1e-4);
            }

            // Case 2: residue 11 (GLU)
            {
                gemmi::Residue* this_res = nullptr;
                for (gemmi::Residue& r : chain.residues) {
                    if (r.seqid.num.value == 11) {
                        this_res = &r;
                        break;
                    }
                }
                ASSERT_NE(this_res, nullptr);

                std::vector<clipper::Coord_orth> residue_atoms = get_atom_coords(*this_res);
                clipper::Coord_orth rc = compute_centroid(residue_atoms);

                // Call our gemmi port
                float r = coot::backrub::residue_radius_gemmi(residue_atoms, rc);
                EXPECT_NEAR(r, 1.54898f, 1e-4);
            }

            break;
        }
    }
    EXPECT_TRUE(found_chain);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}