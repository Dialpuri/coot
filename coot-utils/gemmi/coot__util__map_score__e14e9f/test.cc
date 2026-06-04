#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mtz.hpp>
#include "function.hh"
#include "api/molecules-container.hh"

// Helper to collect non-TER atoms from a gemmi residue
static std::vector<const gemmi::Atom*> collect_atoms(const gemmi::Residue& res) {
    std::vector<const gemmi::Atom*> atoms;
    for (const gemmi::Atom& a : res.atoms) {
        // In gemmi, atoms within a Residue are real atoms (no TER records)
        atoms.push_back(&a);
    }
    return atoms;
}

// Helper to find a residue by chain and seqnum in gemmi
static const gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seqnum) {
    for (gemmi::Chain& ch : model.chains) {
        if (ch.name == chain_id) {
            for (gemmi::Residue& res : ch.residues) {
                if (res.seqid.num.value == seqnum) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, map_score) {
    molecules_container_t mc;
    mc.geometry_init_standard();

    int imol = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(imol, 0);

    int imol_map = mc.read_mtz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
                                "FWT", "PHWT", "", false, false);
    ASSERT_GE(imol_map, 0);

    auto st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    auto& model = st.models[0];

    // Case 1: valid residue (MET at //A/1) — typical case
    {
        const gemmi::Residue* res = find_residue(model, "A", 1);
        ASSERT_NE(res, nullptr);
        EXPECT_STREQ(res->name.c_str(), "MET");
        auto atoms = collect_atoms(*res);
        EXPECT_EQ(atoms.size(), 8u);
        float score = coot::util::map_score_gemmi(atoms, mc[imol_map].xmap, 1);
        EXPECT_NEAR(score, 2.59989, 1e-4);
    }

    // Case 2: another valid residue (PRO at //A/100) — checks different density
    {
        const gemmi::Residue* res = find_residue(model, "A", 100);
        ASSERT_NE(res, nullptr);
        EXPECT_STREQ(res->name.c_str(), "PRO");
        auto atoms = collect_atoms(*res);
        EXPECT_EQ(atoms.size(), 7u);
        float score = coot::util::map_score_gemmi(atoms, mc[imol_map].xmap, 1);
        EXPECT_NEAR(score, 6.88988, 1e-4);
    }

    // Case 3: residue with 0 atomic weighting — tests the weighting parameter path
    {
        const gemmi::Residue* res = find_residue(model, "A", 1);
        ASSERT_NE(res, nullptr);
        EXPECT_STREQ(res->name.c_str(), "MET");
        auto atoms = collect_atoms(*res);
        EXPECT_EQ(atoms.size(), 8u);
        float score = coot::util::map_score_gemmi(atoms, mc[imol_map].xmap, 0);
        EXPECT_NEAR(score, 2.59989, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}