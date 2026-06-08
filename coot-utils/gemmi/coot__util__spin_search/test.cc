#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "api/molecules-container.hh"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, spin_search) {
    {
        // Normal case: ASN 3 with N, CA, CB, CG
        molecules_container_t mc;
        mc.geometry_init_standard();

        int imol_map = mc.read_mtz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
                                   "FWT", "PHWT", "", false, false);
        EXPECT_GE(imol_map, 0);

        clipper::Xmap<float> xmap = mc.get_xmap(imol_map);

        // Load structure with gemmi
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        // Find chain A, residue 3
        gemmi::Chain* chain = nullptr;
        gemmi::Residue* res = nullptr;
        for (gemmi::Chain& ch : st.models[0].chains) {
            if (ch.name == "A") {
                chain = &ch;
                for (gemmi::Residue& r : ch.residues) {
                    if (r.seqid.num.value == 3) {
                        res = &r;
                        break;
                    }
                }
                break;
            }
        }
        EXPECT_NE(res, nullptr);

        // Build torsion with atom specs (names are unpadded in gemmi)
        coot::atom_spec_t atom1("A", 3, "", " N  ", "");
        coot::atom_spec_t atom2("A", 3, "", " CA ", "");
        coot::atom_spec_t atom3("A", 3, "", " CB ", "");
        coot::atom_spec_t atom4("A", 3, "", " CG ", "");

        coot::torsion tors(0, atom1, atom2, atom3, atom4);

        std::pair<float, float> result = coot::util::spin_search_gemmi(xmap, chain, res, tors);

        EXPECT_NEAR(result.first, 24.0f, 1e-4);
        EXPECT_NEAR(result.second, -46.5535f, 1e-4);
        EXPECT_TRUE(result.first > -1000);
    }

    {
        // Edge case: residue without proper atoms (residue 1 lacks CB/CG for this torsion)
        molecules_container_t mc;
        mc.geometry_init_standard();

        int imol_map = mc.read_mtz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
                                   "FWT", "PHWT", "", false, false);
        EXPECT_GE(imol_map, 0);

        clipper::Xmap<float> xmap = mc.get_xmap(imol_map);

        // Load structure with gemmi
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        // Find chain A, residue 1
        gemmi::Chain* chain = nullptr;
        gemmi::Residue* gly = nullptr;
        for (gemmi::Chain& ch : st.models[0].chains) {
            if (ch.name == "A") {
                chain = &ch;
                for (gemmi::Residue& r : ch.residues) {
                    if (r.seqid.num.value == 1) {
                        gly = &r;
                        break;
                    }
                }
                break;
            }
        }
        EXPECT_NE(gly, nullptr);

        coot::atom_spec_t atom1("A", 3, "", "N", "");
        coot::atom_spec_t atom2("A", 3, "", "CA", "");
        coot::atom_spec_t atom3("A", 3, "", "CB", "");
        coot::atom_spec_t atom4("A", 3, "", "CG", "");

        coot::torsion tors_bad(0, atom1, atom2, atom3, atom4);
        std::pair<float, float> result_bad = coot::util::spin_search_gemmi(xmap, chain, gly, tors_bad);

        EXPECT_NEAR(result_bad.first, -1111.1f, 1e-4);
        EXPECT_FALSE(result_bad.first > -1000);
    }
}
