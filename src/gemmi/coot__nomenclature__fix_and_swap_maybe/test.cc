#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "geometry/protein-geometry.hh"

TEST(OracleTest, fix_and_swap_maybe) {
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: apply_swaps = true (actual fixes)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        std::vector<gemmi::CRA> vr = coot::nomenclature::fix_and_swap_maybe_gemmi(st, &geom, true);

        EXPECT_EQ(vr.size(), 8u);

        // residue_0: A 4 PHE
        EXPECT_EQ(std::string(vr[0].chain->name), "A");
        EXPECT_EQ(vr[0].residue->seqid.num.value, 4);
        EXPECT_EQ(std::string(vr[0].residue->name), "PHE");

        // residue_1: A 8 GLU
        EXPECT_EQ(std::string(vr[1].chain->name), "A");
        EXPECT_EQ(vr[1].residue->seqid.num.value, 8);
        EXPECT_EQ(std::string(vr[1].residue->name), "GLU");

        // residue_2: A 19 TYR
        EXPECT_EQ(std::string(vr[2].chain->name), "A");
        EXPECT_EQ(vr[2].residue->seqid.num.value, 19);
        EXPECT_EQ(std::string(vr[2].residue->name), "TYR");

        // residue_3: A 28 GLU
        EXPECT_EQ(std::string(vr[3].chain->name), "A");
        EXPECT_EQ(vr[3].residue->seqid.num.value, 28);
        EXPECT_EQ(std::string(vr[3].residue->name), "GLU");

        // residue_4: A 82 PHE
        EXPECT_EQ(std::string(vr[4].chain->name), "A");
        EXPECT_EQ(vr[4].residue->seqid.num.value, 82);
        EXPECT_EQ(std::string(vr[4].residue->name), "PHE");

        // residue_5: A 208 GLU
        EXPECT_EQ(std::string(vr[5].chain->name), "A");
        EXPECT_EQ(vr[5].residue->seqid.num.value, 208);
        EXPECT_EQ(std::string(vr[5].residue->name), "GLU");

        // residue_6: A 236 TYR
        EXPECT_EQ(std::string(vr[6].chain->name), "A");
        EXPECT_EQ(vr[6].residue->seqid.num.value, 236);
        EXPECT_EQ(std::string(vr[6].residue->name), "TYR");

        // residue_7: A 247 ASP
        EXPECT_EQ(std::string(vr[7].chain->name), "A");
        EXPECT_EQ(vr[7].residue->seqid.num.value, 247);
        EXPECT_EQ(std::string(vr[7].residue->name), "ASP");
    }

    // Case 2: apply_swaps = false (just list, no modifications)
    {
        gemmi::Structure st2 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        std::vector<gemmi::CRA> vr = coot::nomenclature::fix_and_swap_maybe_gemmi(st2, &geom, false);

        EXPECT_EQ(vr.size(), 8u);

        // residue_0: A 4 PHE
        EXPECT_EQ(std::string(vr[0].chain->name), "A");
        EXPECT_EQ(vr[0].residue->seqid.num.value, 4);
        EXPECT_EQ(std::string(vr[0].residue->name), "PHE");

        // residue_1: A 8 GLU
        EXPECT_EQ(std::string(vr[1].chain->name), "A");
        EXPECT_EQ(vr[1].residue->seqid.num.value, 8);
        EXPECT_EQ(std::string(vr[1].residue->name), "GLU");

        // residue_2: A 19 TYR
        EXPECT_EQ(std::string(vr[2].chain->name), "A");
        EXPECT_EQ(vr[2].residue->seqid.num.value, 19);
        EXPECT_EQ(std::string(vr[2].residue->name), "TYR");

        // residue_3: A 28 GLU
        EXPECT_EQ(std::string(vr[3].chain->name), "A");
        EXPECT_EQ(vr[3].residue->seqid.num.value, 28);
        EXPECT_EQ(std::string(vr[3].residue->name), "GLU");

        // residue_4: A 82 PHE
        EXPECT_EQ(std::string(vr[4].chain->name), "A");
        EXPECT_EQ(vr[4].residue->seqid.num.value, 82);
        EXPECT_EQ(std::string(vr[4].residue->name), "PHE");

        // residue_5: A 208 GLU
        EXPECT_EQ(std::string(vr[5].chain->name), "A");
        EXPECT_EQ(vr[5].residue->seqid.num.value, 208);
        EXPECT_EQ(std::string(vr[5].residue->name), "GLU");

        // residue_6: A 236 TYR
        EXPECT_EQ(std::string(vr[6].chain->name), "A");
        EXPECT_EQ(vr[6].residue->seqid.num.value, 236);
        EXPECT_EQ(std::string(vr[6].residue->name), "TYR");

        // residue_7: A 247 ASP
        EXPECT_EQ(std::string(vr[7].chain->name), "A");
        EXPECT_EQ(vr[7].residue->seqid.num.value, 247);
        EXPECT_EQ(std::string(vr[7].residue->name), "ASP");
    }

    // Case 3: empty molecule (edge case — no residues to fix)
    {
        gemmi::Structure st3;

        std::vector<gemmi::CRA> vr = coot::nomenclature::fix_and_swap_maybe_gemmi(st3, &geom, true);

        EXPECT_EQ(vr.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}