#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <clipper/core/xmap.h>
#include <clipper/core/spacegroup.h>
#include <clipper/core/cell.h>
#include <clipper/core/coords.h>
#include "function.hh"

// Helper: find a gemmi::Atom* by CID-like string "//chain/resseq/atomname"
static gemmi::Atom* find_atom_by_cid(const gemmi::Model &model, const std::string &cid) {
    // Parse "//A/10/CA" -> chain="A", resseq=10, atomname="CA"
    size_t pos = 2; // skip "//"
    size_t slash1 = cid.find('/', pos);
    if (slash1 == std::string::npos) return nullptr;
    std::string chain_id = cid.substr(pos, slash1 - pos);

    size_t slash2 = cid.find('/', slash1 + 1);
    if (slash2 == std::string::npos) return nullptr;
    std::string resseq_str = cid.substr(slash1 + 1, slash2 - slash1 - 1);
    int resseq = std::stoi(resseq_str);

    std::string atom_name = cid.substr(slash2 + 1);

    for (const gemmi::Chain &chain : model.chains) {
        if (chain.name != chain_id) continue;
        for (const gemmi::Residue &res : chain.residues) {
            if (res.seqid.num.value != resseq) continue;
            for (const gemmi::Atom &atom : res.atoms) {
                if (atom.name == atom_name) {
                    return const_cast<gemmi::Atom*>(&atom);
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, debug_z_weighted_density_score_new) {
    // Load PDB for atoms
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_FALSE(st.models.empty());
    gemmi::Model &model = st.models[0];

    // Create a minimal valid Xmap for testing density_at_point
    float cell_a = 40.0f, cell_b = 40.0f, cell_c = 80.0f;
    clipper::Cell cell(clipper::Cell_descr(cell_a, cell_b, cell_c, 90, 90, 90));
    clipper::Spacegroup spgr(clipper::Spacegroup::P1);
    clipper::Grid_sampling gs(64, 64, 64);
    clipper::Xmap<float> xmap(spgr, cell, gs);

    // Case 1: A few CA atoms with unit weights
    {
        std::vector<std::pair<gemmi::Atom *, float>> atom_pairs;

        std::string ca10_cid = "//A/10/CA";
        std::string ca20_cid = "//A/20/CA";
        std::string ca30_cid = "//A/30/CA";

        gemmi::Atom *a10 = find_atom_by_cid(model, ca10_cid);
        gemmi::Atom *a20 = find_atom_by_cid(model, ca20_cid);
        gemmi::Atom *a30 = find_atom_by_cid(model, ca30_cid);

        if (a10) { atom_pairs.push_back({a10, 1.0f}); }
        if (a20) { atom_pairs.push_back({a20, 1.0f}); }
        if (a30) { atom_pairs.push_back({a30, 1.0f}); }

        EXPECT_EQ(atom_pairs.size(), 3u);

        // Verify atom properties (gemmi stores unpadded names)
        if (a10) {
            EXPECT_EQ(a10->name, "CA");
            EXPECT_NEAR(a10->pos.x, 27.175, 1e-3);
            EXPECT_NEAR(a10->pos.y, 12.237, 1e-3);
            EXPECT_NEAR(a10->pos.z, 67.71, 1e-2);
        }
        if (a20) {
            EXPECT_EQ(a20->name, "CA");
            EXPECT_NEAR(a20->pos.x, 25.933, 1e-3);
            EXPECT_NEAR(a20->pos.y, 7.675, 1e-3);
            EXPECT_NEAR(a20->pos.z, 71.692, 1e-3);
        }
        if (a30) {
            EXPECT_EQ(a30->name, "CA");
            EXPECT_NEAR(a30->pos.x, 25.226, 1e-3);
            EXPECT_NEAR(a30->pos.y, 1.753, 1e-3);
            EXPECT_NEAR(a30->pos.z, 71.133, 1e-3);
        }

        EXPECT_NO_THROW(coot::util::debug_z_weighted_density_score_new_gemmi(atom_pairs, xmap));
    }

    // Case 2: empty vector (edge case)
    {
        std::vector<std::pair<gemmi::Atom *, float>> atom_pairs;

        EXPECT_EQ(atom_pairs.size(), 0u);

        EXPECT_NO_THROW(coot::util::debug_z_weighted_density_score_new_gemmi(atom_pairs, xmap));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
