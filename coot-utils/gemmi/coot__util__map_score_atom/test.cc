#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "coot-utils/coot-map-utils.hh"

// Helper: find atom by CID string like "//A/50/CA" in a gemmi::Structure
static gemmi::Atom* find_atom_by_cid(const std::string &cid, gemmi::Structure &st) {
    // Parse CID: "//<chain>/<res_seq>/<atom_name>"
    // Skip leading "//"
    std::string rest = cid.substr(2);
    auto p1 = rest.find('/');
    if (p1 == std::string::npos) return nullptr;
    std::string chain_id = rest.substr(0, p1);
    rest = rest.substr(p1 + 1);
    auto p2 = rest.find('/');
    if (p2 == std::string::npos) return nullptr;
    std::string res_seq_str = rest.substr(0, p2);
    std::string atom_name = rest.substr(p2 + 1);

    int res_seq = std::stoi(res_seq_str);

    gemmi::Model &model = st.models[0];
    for (gemmi::Chain &chain : model.chains) {
        if (chain.name != chain_id) continue;
        for (gemmi::Residue &res : chain.residues) {
            if (res.seqid.num.value != res_seq) continue;
            for (gemmi::Atom &atom : res.atoms) {
                // Trim atom name for comparison
                std::string trimmed = atom.name;
                trimmed.erase(0, trimmed.find_first_not_of(" "));
                trimmed.erase(trimmed.find_last_not_of(" ") + 1);
                if (trimmed == atom_name) return &atom;
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, map_score_atom) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    clipper::Xmap<float> xmap;
    bool r = coot::util::map_fill_from_mtz(&xmap,
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
        "FWT", "PHWT", "", false, 1.8);
    EXPECT_TRUE(r);

    // --- Case 1: valid atom at residue A/50 CA ---
    {
        std::string cid = "//A/50/CA";
        gemmi::Atom *atom = find_atom_by_cid(cid, st);
        EXPECT_NE(atom, nullptr);
        EXPECT_NEAR(atom->pos.x, 8.709, 1e-4);
        EXPECT_NEAR(atom->pos.y, -3.609, 1e-4);
        EXPECT_NEAR(atom->pos.z, 56.949, 1e-4);
        float score = coot::util::map_score_atom_gemmi(atom, xmap);
        EXPECT_NEAR(score, 0.412424, 1e-4);
    }

    // --- Case 2: null atom pointer (edge case) ---
    {
        std::string cid = "//Z/9999/CA";
        gemmi::Atom *atom = find_atom_by_cid(cid, st);
        EXPECT_EQ(atom, nullptr);
        float score = coot::util::map_score_atom_gemmi(atom, xmap);
        EXPECT_NEAR(score, 0.0, 1e-4);
    }

    // --- Case 3: another valid atom ---
    {
        std::string cid = "//A/1/CA";
        gemmi::Atom *atom = find_atom_by_cid(cid, st);
        EXPECT_NE(atom, nullptr);
        float score = coot::util::map_score_atom_gemmi(atom, xmap);
        EXPECT_NEAR(score, 0.566512, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}