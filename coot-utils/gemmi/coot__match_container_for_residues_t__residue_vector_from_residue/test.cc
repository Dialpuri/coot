#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <algorithm>

TEST(OracleTest, residue_vector_from_residue) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // Find chain A
    gemmi::Chain* chainA = nullptr;
    for (auto& ch : model.chains) {
        if (ch.name == "A") {
            chainA = &ch;
            break;
        }
    }
    ASSERT_NE(chainA, nullptr);

    coot::match_container_for_residues_t container;

    // Case 1: residue in the middle of chain A (seq 50) — should find contiguous neighbours
    {
        gemmi::Residue* res = nullptr;
        for (auto& r : chainA->residues) {
            if (r.seqid.num.value == 50) {
                res = &r;
                break;
            }
        }
        // fallback to first real residue
        if (!res) {
            res = &chainA->residues[0];
        }
        gemmi::CRA cra{chainA, res, nullptr};
        std::vector<gemmi::Residue*> v = container.residue_vector_from_residue_gemmi(cra);
        EXPECT_EQ(v.size(), 27u);
        if (!v.empty()) {
            std::sort(v.begin(), v.end(), [](gemmi::Residue* a, gemmi::Residue* b) {
                return a->seqid.num.value < b->seqid.num.value;
            });
            EXPECT_EQ(v.front()->seqid.num.value, 46);
            EXPECT_EQ(v.back()->seqid.num.value, 72);
        }
    }

    // Case 2: first residue in chain A (seq 1) — should still find contiguous but only forward
    {
        gemmi::Residue* res = &chainA->residues[0];
        gemmi::CRA cra{chainA, res, nullptr};
        std::vector<gemmi::Residue*> v = container.residue_vector_from_residue_gemmi(cra);
        EXPECT_EQ(v.size(), 36u);
        if (!v.empty()) {
            std::sort(v.begin(), v.end(), [](gemmi::Residue* a, gemmi::Residue* b) {
                return a->seqid.num.value < b->seqid.num.value;
            });
            EXPECT_EQ(v.front()->seqid.num.value, 1);
            EXPECT_EQ(v.back()->seqid.num.value, 36);
        }
    }

    // Case 3: last residue in chain A — should find contiguous backward
    {
        gemmi::Residue* res = &chainA->residues.back();
        gemmi::CRA cra{chainA, res, nullptr};
        std::vector<gemmi::Residue*> v = container.residue_vector_from_residue_gemmi(cra);
        EXPECT_EQ(v.size(), 133u);
        if (!v.empty()) {
            std::sort(v.begin(), v.end(), [](gemmi::Residue* a, gemmi::Residue* b) {
                return a->seqid.num.value < b->seqid.num.value;
            });
            EXPECT_EQ(v.front()->seqid.num.value, 166);
            EXPECT_EQ(v.back()->seqid.num.value, 298);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}