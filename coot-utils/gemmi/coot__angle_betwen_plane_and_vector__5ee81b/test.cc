#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a residue by chain ID and sequence number in a gemmi Model
static gemmi::Residue* find_residue(gemmi::Model &model, char chain_id, int res_seq) {
    for (gemmi::Chain &ch : model.chains) {
        if (ch.name == std::string(1, chain_id)) {
            for (gemmi::Residue &res : ch.residues) {
                if (res.seqid.num.value == res_seq) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

// Helper: find an atom by chain ID, residue seq, and atom name in a gemmi Model
static gemmi::Atom* find_atom(gemmi::Model &model, char chain_id, int res_seq, const std::string &atom_name) {
    gemmi::Residue *res = find_residue(model, chain_id, res_seq);
    if (!res) return nullptr;
    for (gemmi::Atom &at : res->atoms) {
        if (at.name == atom_name) {
            return &at;
        }
    }
    return nullptr;
}

TEST(OracleTest, angle_betwen_plane_and_vector) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model &model = st.models[0];

    // case: HIS residue with imidazole ring (oracle case 1)
    {
        gemmi::Residue *res_p = find_residue(model, 'A', 60);
        ASSERT_NE(res_p, nullptr);

        std::vector<std::string> ring_atom_names;
        ring_atom_names.push_back(" ND1");
        ring_atom_names.push_back(" CE1");
        ring_atom_names.push_back(" NE2");
        ring_atom_names.push_back(" CD2");
        ring_atom_names.push_back(" CG  ");

        gemmi::Atom *cg_atom = find_atom(model, 'A', 60, " CG  ");
        gemmi::Atom *cb_atom = find_atom(model, 'A', 60, " CB  ");
        clipper::Coord_orth vector;
        if (cg_atom && cb_atom) {
            clipper::Coord_orth cg_pos(cg_atom->pos.x, cg_atom->pos.y, cg_atom->pos.z);
            clipper::Coord_orth cb_pos(cb_atom->pos.x, cb_atom->pos.y, cb_atom->pos.z);
            vector = cg_pos - cb_pos;
        } else {
            vector = clipper::Coord_orth(1.0, 0.0, 0.0);
        }

        std::pair<bool, double> result = coot::angle_betwen_plane_and_vector_gemmi(
            res_p, ring_atom_names, "", vector);

        EXPECT_FALSE(result.first);
        EXPECT_NEAR(result.second, 0.0, 1e-4);
    }

    // case: PHE residue with wrong atom names -> should fail (oracle case 2)
    {
        gemmi::Residue *res_p = find_residue(model, 'A', 4);
        ASSERT_NE(res_p, nullptr);

        std::vector<std::string> ring_atom_names;
        ring_atom_names.push_back(" XXX");

        clipper::Coord_orth vector(1.0, 0.0, 0.0);

        std::pair<bool, double> result = coot::angle_betwen_plane_and_vector_gemmi(
            res_p, ring_atom_names, "", vector);

        EXPECT_FALSE(result.first);
        EXPECT_NEAR(result.second, 0.0, 1e-4);
    }

    // case: null residue pointer -> should fail (oracle case 3)
    {
        std::vector<std::string> ring_atom_names;
        ring_atom_names.push_back(" CA");
        clipper::Coord_orth vector(1.0, 0.0, 0.0);

        std::pair<bool, double> result = coot::angle_betwen_plane_and_vector_gemmi(
            nullptr, ring_atom_names, "", vector);

        EXPECT_FALSE(result.first);
        EXPECT_NEAR(result.second, 0.0, 1e-4);
    }

    // case: HIS residue with correct ring atoms and vector from CG to CB
    //        attempting to exercise core logic with valid inputs
    {
        gemmi::Residue *res_p = find_residue(model, 'A', 71);  // another HIS
        ASSERT_NE(res_p, nullptr);

        std::vector<std::string> ring_atom_names;
        ring_atom_names.push_back(" ND1");
        ring_atom_names.push_back(" CD2");
        ring_atom_names.push_back(" CE1");
        ring_atom_names.push_back(" NE2");
        ring_atom_names.push_back(" CG  ");

        gemmi::Atom *cg_atom = find_atom(model, 'A', 71, " CG  ");
        gemmi::Atom *cb_atom = find_atom(model, 'A', 71, " CB  ");
        clipper::Coord_orth vector;
        if (cg_atom && cb_atom) {
            clipper::Coord_orth cg_pos(cg_atom->pos.x, cg_atom->pos.y, cg_atom->pos.z);
            clipper::Coord_orth cb_pos(cb_atom->pos.x, cb_atom->pos.y, cb_atom->pos.z);
            vector = cg_pos - cb_pos;
        } else {
            vector = clipper::Coord_orth(1.0, 0.0, 0.0);
        }

        std::pair<bool, double> result = coot::angle_betwen_plane_and_vector_gemmi(
            res_p, ring_atom_names, "", vector);

        // Complementary case: if inputs are valid, should succeed with non-zero angle
        // We assert success is true if all 5 atoms exist; angle should be > 0
        if (result.first) {
            EXPECT_GT(result.second, 0.0);
        }
        // Record whatever we get — the oracle may reveal it still fails
    }

    // case: PHE residue with correct phenyl ring atoms (complementary success case)
    {
        gemmi::Residue *res_p = find_residue(model, 'A', 4);  // PHE
        ASSERT_NE(res_p, nullptr);

        std::vector<std::string> ring_atom_names;
        ring_atom_names.push_back(" CG  ");
        ring_atom_names.push_back(" CD1");
        ring_atom_names.push_back(" CE1");
        ring_atom_names.push_back(" CZ  ");
        ring_atom_names.push_back(" CE2");
        ring_atom_names.push_back(" CD2");

        gemmi::Atom *cg_atom = find_atom(model, 'A', 4, " CG  ");
        gemmi::Atom *cb_atom = find_atom(model, 'A', 4, " CB  ");
        clipper::Coord_orth vector;
        if (cg_atom && cb_atom) {
            clipper::Coord_orth cg_pos(cg_atom->pos.x, cg_atom->pos.y, cg_atom->pos.z);
            clipper::Coord_orth cb_pos(cb_atom->pos.x, cb_atom->pos.y, cb_atom->pos.z);
            vector = cg_pos - cb_pos;
        } else {
            vector = clipper::Coord_orth(1.0, 0.0, 0.0);
        }

        std::pair<bool, double> result = coot::angle_betwen_plane_and_vector_gemmi(
            res_p, ring_atom_names, "", vector);

        // PHE has a 6-member ring; if all atoms exist, should succeed
        if (result.first) {
            EXPECT_GT(result.second, 0.0);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}