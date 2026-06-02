#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static std::string pad4(const std::string &s) {
    std::string r(s);
    while (r.size() < 4) r = " " + r;
    return r;
}

TEST(OracleTest, test_and_fix_ASP_GLU_nomenclature_errors) {
    // Setup: load PDB and find residues
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model &model = st.models[0];
    gemmi::Chain &chain = model.chains[0];

    // Find an ASP residue and a GLU residue
    gemmi::Residue *asp_res = nullptr;
    gemmi::Residue *glu_res = nullptr;

    for (gemmi::Residue &res : chain.residues) {
        std::string rname = res.name;
        if (rname == "ASP" && !asp_res) {
            asp_res = &res;
        }
        if (rname == "GLU" && !glu_res) {
            glu_res = &res;
        }
        if (asp_res && glu_res) break;
    }

    // Case 1: ASP residue (apply_swap_if_found = true)
    {
        if (asp_res) {
            // Get atom names before (pad to 4 chars to match MMDB convention)
            std::string od1_before, od2_before;
            for (gemmi::Atom &atom : asp_res->atoms) {
                if (atom.name == "OD1") od1_before = pad4(atom.name);
                if (atom.name == "OD2") od2_before = pad4(atom.name);
            }

            // Oracle ground truth
            EXPECT_EQ(od1_before, " OD1");
            EXPECT_EQ(od2_before, " OD2");

            int result = coot::nomenclature::test_and_fix_ASP_GLU_nomenclature_errors_gemmi(*asp_res, true);

            // Get atom names after (pad to 4 chars to match MMDB convention)
            std::string od1_after, od2_after;
            for (gemmi::Atom &atom : asp_res->atoms) {
                if (atom.name == "OD1") od1_after = pad4(atom.name);
                if (atom.name == "OD2") od2_after = pad4(atom.name);
            }

            // Oracle ground truth
            EXPECT_EQ(result, 0);
            EXPECT_EQ(od1_after, " OD1");
            EXPECT_EQ(od2_after, " OD2");
        }
    }

    // Case 2: GLU residue (apply_swap_if_found = false — just detect)
    {
        if (glu_res) {
            int result = coot::nomenclature::test_and_fix_ASP_GLU_nomenclature_errors_gemmi(*glu_res, false);
            // Oracle ground truth
            EXPECT_EQ(result, 0);
        }
    }

    // Case 3: Non-ASP/GLU residue (should return 0, no-op)
    {
        if (!chain.residues.empty()) {
            gemmi::Residue &other_res = chain.residues[0];
            std::string other_name = other_res.name;
            if (other_name != "ASP" && other_name != "GLU") {
                int result = coot::nomenclature::test_and_fix_ASP_GLU_nomenclature_errors_gemmi(other_res, true);
                // Oracle ground truth
                EXPECT_EQ(result, 0);
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}