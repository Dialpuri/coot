#include <gtest/gtest.h>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, loop_path) {
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");

    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    gemmi::Chain* chain = nullptr;
    for (auto& ch : model.chains) {
        if (ch.name == "A") { chain = &ch; break; }
    }
    ASSERT_NE(chain, nullptr);

    // --- Case 1: Valid 4 atoms (consecutive residues 1-4, CA atoms) ---
    {
        ASSERT_GE(chain->residues.size(), 4u);
        gemmi::Residue& r1 = chain->residues[0];   // seqnum 1
        gemmi::Residue& r2 = chain->residues[1];   // seqnum 2
        gemmi::Residue& r3 = chain->residues[2];   // seqnum 3
        gemmi::Residue& r4 = chain->residues[3];   // seqnum 4

        // Find CA atoms in each residue
        const gemmi::Atom* start_back_2 = nullptr;
        const gemmi::Atom* start        = nullptr;
        const gemmi::Atom* end          = nullptr;
        const gemmi::Atom* end_plus_2   = nullptr;

        for (const gemmi::Atom& atom : r1.atoms) {
            if (atom.name == "CA") { start_back_2 = &atom; break; }
        }
        for (const gemmi::Atom& atom : r2.atoms) {
            if (atom.name == "CA") { start = &atom; break; }
        }
        for (const gemmi::Atom& atom : r3.atoms) {
            if (atom.name == "CA") { end = &atom; break; }
        }
        for (const gemmi::Atom& atom : r4.atoms) {
            if (atom.name == "CA") { end_plus_2 = &atom; break; }
        }

        gemmi::CRA cra_s_back_2 = {chain, &r1, const_cast<gemmi::Atom*>(start_back_2)};
        gemmi::CRA cra_s        = {chain, &r2, const_cast<gemmi::Atom*>(start)};
        gemmi::CRA cra_e        = {chain, &r3, const_cast<gemmi::Atom*>(end)};
        gemmi::CRA cra_e_plus_2 = {chain, &r4, const_cast<gemmi::Atom*>(end_plus_2)};

        unsigned int n_line_segments = 8;

        auto result = coot::loop_path_gemmi(cra_s_back_2, cra_s, cra_e, cra_e_plus_2, n_line_segments);

        EXPECT_TRUE(result.first);
        EXPECT_EQ(result.second.size(), 16u);

        if (!result.second.empty()) {
            const coot::CartesianPair &first = result.second.front();
            const coot::CartesianPair &last = result.second.back();

            EXPECT_NEAR(first.getStart().x(), 13.627, 1e-4);
            EXPECT_NEAR(first.getStart().y(), 4.716, 1e-4);
            EXPECT_NEAR(first.getStart().z(), 78.375, 1e-4);

            EXPECT_NEAR(last.getFinish().x(), 15.7583, 1e-4);
            EXPECT_NEAR(last.getFinish().y(), 1.75744, 1e-4);
            EXPECT_NEAR(last.getFinish().z(), 78.3895, 1e-4);
        }
    }

    // --- Case 2: Null pointer guard (early return with blank) ---
    {
        gemmi::CRA empty = {nullptr, nullptr, nullptr};
        auto result = coot::loop_path_gemmi(empty, empty, empty, empty, 8);

        EXPECT_FALSE(result.first);
        EXPECT_EQ(result.second.size(), 0u);
    }
}