#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cstdlib>

#include "geometry/protein-geometry.hh"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTestGemmi, make_rama_plot_restraints_ng) {
    // Shared setup — load PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());
    gemmi::Chain& chain = model.chains[0];
    ASSERT_FALSE(chain.residues.empty());

    coot::protein_geometry geom;
    geom.init_standard();

    // --- Case 1: Original oracle case (consecutive standard amino acids) ---
    {
        // Collect indices of standard amino acid residues that are
        // consecutive in the chain
        std::vector<int> std_aa_indices;
        for (int i = 0; i < static_cast<int>(chain.residues.size()); i++) {
            std::string rn = chain.residues[i].name;
            if (coot::util::is_standard_amino_acid_name(rn)) {
                std_aa_indices.push_back(i);
            }
        }

        // Find 3 consecutive indices
        std::vector<int> consecutive_indices;
        for (size_t i = 1; i + 1 < std_aa_indices.size(); i++) {
            if ((std_aa_indices[i] - std_aa_indices[i - 1]) == 1 &&
                (std_aa_indices[i + 1] - std_aa_indices[i]) == 1) {
                consecutive_indices.push_back(std_aa_indices[i - 1]);
                consecutive_indices.push_back(std_aa_indices[i]);
                consecutive_indices.push_back(std_aa_indices[i + 1]);
                break;
            }
        }

        EXPECT_GE(consecutive_indices.size(), 3u);

        // Verify oracle input values (seq nums only)
        EXPECT_EQ(chain.residues[consecutive_indices[0]].seqid.num.value, 1);
        EXPECT_EQ(chain.residues[consecutive_indices[1]].seqid.num.value, 2);
        EXPECT_EQ(chain.residues[consecutive_indices[2]].seqid.num.value, 3);

        // Get pointers for map/set building
        gemmi::Residue* r0 = &chain.residues[consecutive_indices[0]];
        gemmi::Residue* r1 = &chain.residues[consecutive_indices[1]];
        gemmi::Residue* r2 = &chain.residues[consecutive_indices[2]];

        // Build residues list as CRA (carries chain context)
        std::vector<gemmi::CRA> residues;
        std::vector<bool> fixed_flags;
        residues.push_back(gemmi::CRA{&chain, r0, nullptr});
        residues.push_back(gemmi::CRA{&chain, r1, nullptr});
        residues.push_back(gemmi::CRA{&chain, r2, nullptr});
        fixed_flags = {false, false, false};  // all not fixed
        EXPECT_EQ(residues.size(), 3u);

        // Build residue_link_vector_map: for the middle residue, map to its neighbors
        std::map<gemmi::Residue*, std::vector<gemmi::Residue*>> residue_link_vector_map;
        residue_link_vector_map[r1] = {r0, r2};

        // Build residue_pair_link_set
        std::set<std::pair<gemmi::Residue*, gemmi::Residue*>> residue_pair_link_set;
        residue_pair_link_set.insert({r0, r1});
        residue_pair_link_set.insert({r1, r2});

        // Call the ported function
        auto triples = coot::make_rama_plot_restraints_ng_gemmi(
            residues, fixed_flags, residue_link_vector_map, residue_pair_link_set, geom);

        // Should produce exactly 1 rama triple (original: n_after - n_before == 1)
        EXPECT_EQ(triples.size(), 1u);

        if (!triples.empty()) {
            // Verify the triple's link type is TRANS or PTRANS
            EXPECT_TRUE(triples[0].link_type == "TRANS" || triples[0].link_type == "PTRANS");
            // Verify residues are the ones we passed in
            EXPECT_EQ(triples[0].r_1, r0);
            EXPECT_EQ(triples[0].r_2, r1);
            EXPECT_EQ(triples[0].r_3, r2);
        }
    }

    // --- Case 2: Empty residue_link_vector_map should add 0 restraints ---
    {
        // Use first 3 standard residues by index
        std::vector<int> std_aa_indices;
        for (int i = 0; i < static_cast<int>(chain.residues.size()); i++) {
            std::string rn = chain.residues[i].name;
            if (coot::util::is_standard_amino_acid_name(rn)) {
                std_aa_indices.push_back(i);
                if (std_aa_indices.size() >= 3) break;
            }
        }

        EXPECT_GE(std_aa_indices.size(), 3u);

        gemmi::Residue* r0 = &chain.residues[std_aa_indices[0]];
        gemmi::Residue* r1 = &chain.residues[std_aa_indices[1]];
        gemmi::Residue* r2 = &chain.residues[std_aa_indices[2]];

        std::vector<gemmi::CRA> residues;
        std::vector<bool> fixed_flags;
        residues.push_back(gemmi::CRA{&chain, r0, nullptr});
        residues.push_back(gemmi::CRA{&chain, r1, nullptr});
        residues.push_back(gemmi::CRA{&chain, r2, nullptr});
        fixed_flags = {false, false, false};

        // Empty residue_link_vector_map — no residues to make rama restraints for
        std::map<gemmi::Residue*, std::vector<gemmi::Residue*>> residue_link_vector_map;
        std::set<std::pair<gemmi::Residue*, gemmi::Residue*>> residue_pair_link_set;

        auto triples = coot::make_rama_plot_restraints_ng_gemmi(
            residues, fixed_flags, residue_link_vector_map, residue_pair_link_set, geom);

        // Empty map means no restraints added
        EXPECT_EQ(triples.size(), 0u);
    }

    // --- Case 3: Non-standard residue should not add restraints ---
    {
        // Find a residue that is NOT a standard amino acid
        std::vector<int> non_std_indices;
        for (int i = 0; i < static_cast<int>(chain.residues.size()); i++) {
            std::string rn = chain.residues[i].name;
            if (!coot::util::is_standard_amino_acid_name(rn)) {
                non_std_indices.push_back(i);
                if (non_std_indices.size() >= 1) break;
            }
        }

        if (non_std_indices.size() >= 1) {
            gemmi::Residue* ns_res = &chain.residues[non_std_indices[0]];

            // Build residues with a non-standard residue
            std::vector<gemmi::CRA> residues;
            std::vector<bool> fixed_flags;
            residues.push_back(gemmi::CRA{&chain, ns_res, nullptr});
            fixed_flags.push_back(false);

            // Build residue_link_vector_map with non-standard residue, no neighbors
            std::map<gemmi::Residue*, std::vector<gemmi::Residue*>> residue_link_vector_map;
            residue_link_vector_map[ns_res] = {};

            std::set<std::pair<gemmi::Residue*, gemmi::Residue*>> residue_pair_link_set;

            auto triples = coot::make_rama_plot_restraints_ng_gemmi(
                residues, fixed_flags, residue_link_vector_map, residue_pair_link_set, geom);

            // Non-standard residue with no neighbors (and only 1 residue, < 3) should add 0
            EXPECT_EQ(triples.size(), 0u);
        }
    }

    // --- Case 4: All-fixed residues should not add restraints ---
    {
        // Use the same consecutive residues as Case 1
        std::vector<int> std_aa_indices;
        for (int i = 0; i < static_cast<int>(chain.residues.size()); i++) {
            std::string rn = chain.residues[i].name;
            if (coot::util::is_standard_amino_acid_name(rn)) {
                std_aa_indices.push_back(i);
            }
            if (std_aa_indices.size() >= 6) break;
        }

        std::vector<int> consecutive_indices;
        for (size_t i = 1; i + 1 < std_aa_indices.size(); i++) {
            if ((std_aa_indices[i] - std_aa_indices[i - 1]) == 1 &&
                (std_aa_indices[i + 1] - std_aa_indices[i]) == 1) {
                consecutive_indices.push_back(std_aa_indices[i - 1]);
                consecutive_indices.push_back(std_aa_indices[i]);
                consecutive_indices.push_back(std_aa_indices[i + 1]);
                break;
            }
        }

        if (consecutive_indices.size() >= 3) {
            gemmi::Residue* r0 = &chain.residues[consecutive_indices[0]];
            gemmi::Residue* r1 = &chain.residues[consecutive_indices[1]];
            gemmi::Residue* r2 = &chain.residues[consecutive_indices[2]];

            std::vector<gemmi::CRA> residues;
            residues.push_back(gemmi::CRA{&chain, r0, nullptr});
            residues.push_back(gemmi::CRA{&chain, r1, nullptr});
            residues.push_back(gemmi::CRA{&chain, r2, nullptr});
            // ALL FIXED
            std::vector<bool> fixed_flags = {true, true, true};

            std::map<gemmi::Residue*, std::vector<gemmi::Residue*>> residue_link_vector_map;
            residue_link_vector_map[r1] = {r0, r2};

            std::set<std::pair<gemmi::Residue*, gemmi::Residue*>> residue_pair_link_set;
            residue_pair_link_set.insert({r0, r1});
            residue_pair_link_set.insert({r1, r2});

            auto triples = coot::make_rama_plot_restraints_ng_gemmi(
                residues, fixed_flags, residue_link_vector_map, residue_pair_link_set, geom);

            // All fixed => skip, 0 restraints
            EXPECT_EQ(triples.size(), 0u);
        }
    }
}
