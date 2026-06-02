#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static const char* example_pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";

TEST(OracleTest, get_base_normal) {
    gemmi::Structure st = gemmi::read_pdb_file(example_pdb_path);
    gemmi::Model& model = st.models[0];

    // Case 1: A protein residue (ILE at pos 10) — no base atoms, should return false
    {
        gemmi::Residue* res = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& r : chain.residues) {
                    if (r.seqid.num.value == 10) {
                        res = &r;
                        break;
                    }
                }
            }
        }
        ASSERT_NE(res, nullptr);

        std::pair<bool, clipper::Coord_orth> result = coot::stack_and_pair::get_base_normal_gemmi(res);

        EXPECT_FALSE(result.first);
    }

    // Case 2: Another protein residue (MET at pos 1) — also no base atoms, should return false
    {
        gemmi::Residue* res = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& r : chain.residues) {
                    if (r.seqid.num.value == 1) {
                        res = &r;
                        break;
                    }
                }
            }
        }
        ASSERT_NE(res, nullptr);

        std::pair<bool, clipper::Coord_orth> result = coot::stack_and_pair::get_base_normal_gemmi(res);

        EXPECT_FALSE(result.first);
    }

    // Case 3: Nucleic acid residue from chain C — may have base atoms (C4', C1', N9/N1)
    {
        gemmi::Residue* res = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "C") {
                for (gemmi::Residue& r : chain.residues) {
                    if (r.seqid.num.value == 1) {
                        res = &r;
                        break;
                    }
                }
            }
        }
        if (res != nullptr) {
            std::pair<bool, clipper::Coord_orth> result = coot::stack_and_pair::get_base_normal_gemmi(res);
            EXPECT_NO_THROW(result);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}