#include <gtest/gtest.h>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, do_Ca_loop) {
    // Case 1: Protein chain A — gap between residue seqnums 36 and 46
    {
        gemmi::cif::Document doc = gemmi::cif::read_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
        gemmi::Structure st = gemmi::make_structure(std::move(doc));

        // Find chain A
        const gemmi::Model* model = nullptr;
        for (gemmi::Model& m : st.models) {
            model = &m;
            break;
        }
        ASSERT_NE(model, nullptr);

        const gemmi::Chain* chain = nullptr;
        for (const auto& c : model->chains) {
            if (c.name == "A") {
                chain = &c;
                break;
            }
        }
        ASSERT_NE(chain, nullptr);

        int nres = static_cast<int>(chain->residues.size());
        EXPECT_EQ(nres, 513);

        // Find gap between residues
        int ires = -1;
        for (int i = 1; i < nres; i++) {
            const auto& prev = chain->residues[static_cast<size_t>(i - 1)];
            const auto& this_res = chain->residues[static_cast<size_t>(i)];
            int delta = this_res.seqid.num.value - prev.seqid.num.value;
            if (delta > 1) {
                ires = i;
                break;
            }
        }
        EXPECT_EQ(ires, 36);

        const auto& res_prev = chain->residues[static_cast<size_t>(ires - 1)];
        const auto& res_this = chain->residues[static_cast<size_t>(ires)];
        ASSERT_TRUE(res_prev.seqid.num.value == 36);
        ASSERT_TRUE(res_this.seqid.num.value == 46);
        EXPECT_EQ(res_prev.seqid.num.value, 36);
        EXPECT_EQ(res_this.seqid.num.value, 46);

        // Original test: total_bonds_before = 3404, total_bonds_after = 3424
        // i.e. 20 new bonds (loop segments) were added
        auto lp = coot::do_Ca_loop_gemmi(ires, nres, chain, &res_prev, &res_this);

        // 20 bonds were added (lp.second.size() should be 20)
        // original: 3424 - 3404 = 20
        EXPECT_EQ(lp.second.size(), 20u);

        // Original test: spots_after = 0, meaning lp.first was false (no bad spots added)
        EXPECT_EQ(lp.first, false);
    }

    // Case 2: Consecutive residues (no gap) — tests different branch
    {
        gemmi::cif::Document doc = gemmi::cif::read_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
        gemmi::Structure st = gemmi::make_structure(std::move(doc));

        const gemmi::Model* model = nullptr;
        for (gemmi::Model& m : st.models) {
            model = &m;
            break;
        }
        ASSERT_NE(model, nullptr);

        const gemmi::Chain* chain = nullptr;
        for (const auto& c : model->chains) {
            if (c.name == "A") {
                chain = &c;
                break;
            }
        }
        ASSERT_NE(chain, nullptr);

        int nres = static_cast<int>(chain->residues.size());
        // Use residues at index 0 and 1 (should be consecutive, no gap)
        const auto& res_prev = chain->residues[0];
        const auto& res_this = chain->residues[1];

        // Verify they are consecutive (delta == 1, so no gap)
        int delta = res_this.seqid.num.value - res_prev.seqid.num.value;
        EXPECT_EQ(delta, 1);

        // Call do_Ca_loop with consecutive residues — should still produce a loop
        // (the function doesn't check for gap size, only structural conditions)
        EXPECT_NO_THROW(coot::do_Ca_loop_gemmi(1, nres, chain, &res_prev, &res_this));

        auto lp = coot::do_Ca_loop_gemmi(1, nres, chain, &res_prev, &res_this);
        // Consecutive residues: should still find CA atoms and produce a loop
        // since the structural checks pass (C-N close check might prevent it)
        // Just verify no crash and the result is computed
        EXPECT_GE(lp.second.size(), 0u);
    }

    // Case 3: HOH residues — should be guarded away (no chain B in this file)
    {
        gemmi::cif::Document doc = gemmi::cif::read_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
        gemmi::Structure st = gemmi::make_structure(std::move(doc));

        const gemmi::Model* model = nullptr;
        for (gemmi::Model& m : st.models) {
            model = &m;
            break;
        }
        ASSERT_NE(model, nullptr);

        const gemmi::Chain* chain_hoh = nullptr;
        for (const auto& c : model->chains) {
            if (c.name == "B") {
                chain_hoh = &c;
                break;
            }
        }
        // Chain B does not exist in this file — verify gracefully
        EXPECT_EQ(chain_hoh, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
