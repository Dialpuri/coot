#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <algorithm>
#include <cmath>

TEST(OracleTest, typed_distances_generate) {
    // --- Setup: load PDB ---
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Collect all non-hydrogen atoms as CRAs (equivalent to MMDB selection of !H)
    std::vector<gemmi::CRA> selected_atoms;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            for (gemmi::Atom& atom : res.atoms) {
                if (!atom.element.is_hydrogen()) {
                    selected_atoms.push_back({&chain, &res, &atom});
                }
            }
        }
    }
    EXPECT_EQ(selected_atoms.size(), 2107u);

    // Build residues_within_window
    std::map<gemmi::CRA, std::vector<gemmi::CRA>, coot::typed_distances::cra_less> residues_within_window;

    // Build residue_distances_map
    std::map<gemmi::CRA, std::map<int, std::vector<float>>, coot::typed_distances::cra_less> residue_distances_map;

    // First call — simulate auto-init
    coot::typed_distances::generate_gemmi(model, selected_atoms, residue_distances_map, residues_within_window);

    size_t after_init_size = residue_distances_map.size();
    EXPECT_EQ(after_init_size, 268u);

    // Clear the map so we can demonstrate generate() directly
    residue_distances_map.clear();
    size_t before_generate_size = residue_distances_map.size();
    EXPECT_EQ(before_generate_size, 0u);

    // Call generate — this is the function under test
    coot::typed_distances::generate_gemmi(model, selected_atoms, residue_distances_map, residues_within_window);

    size_t after_generate_size = residue_distances_map.size();
    EXPECT_EQ(after_generate_size, 268u);

    // Verify sample residue distances from the map
    // We look up residues by identity (chain, seqnum, resname) rather than
    // relying on iteration order, since gemmi's cra_less (pointer-based)
    // orders differently from MMDB's pointer allocation order.

    // Build a lookup from (chain_id, seqnum) -> CRA for chain A
    struct res_key {
        std::string chain_id;
        int seqnum;
        bool operator<(const res_key& o) const {
            if (chain_id != o.chain_id) return chain_id < o.chain_id;
            return seqnum < o.seqnum;
        }
    };
    std::map<res_key, gemmi::CRA> res_lookup;
    for (auto const& outer : residue_distances_map) {
        gemmi::CRA cra = outer.first;
        ASSERT_NE(cra.chain, nullptr);
        ASSERT_NE(cra.residue, nullptr);
        res_key rk{cra.chain->name, cra.residue->seqid.num.value};
        res_lookup[rk] = cra;
    }

    // Sample 1: chain A, seqnum 6, LYS -> 102 total distances
    {
        auto it = res_lookup.find(res_key{"A", 6});
        ASSERT_NE(it, res_lookup.end());
        gemmi::CRA cra = it->second;
        EXPECT_EQ(cra.chain->name, "A");
        EXPECT_EQ(cra.residue->seqid.num.value, 6);
        EXPECT_EQ(cra.residue->name, "LYS");
        int total_dists = 0;
        for (auto const& outer : residue_distances_map) {
            if (outer.first.residue == cra.residue) {
                for (auto const& inner : outer.second) {
                    total_dists += static_cast<int>(inner.second.size());
                }
            }
        }
        EXPECT_EQ(total_dists, 102);
    }

    // Sample 2: chain A, seqnum 5, GLN -> 114 total distances
    {
        auto it = res_lookup.find(res_key{"A", 5});
        ASSERT_NE(it, res_lookup.end());
        gemmi::CRA cra = it->second;
        EXPECT_EQ(cra.chain->name, "A");
        EXPECT_EQ(cra.residue->seqid.num.value, 5);
        EXPECT_EQ(cra.residue->name, "GLN");
        int total_dists = 0;
        for (auto const& outer : residue_distances_map) {
            if (outer.first.residue == cra.residue) {
                for (auto const& inner : outer.second) {
                    total_dists += static_cast<int>(inner.second.size());
                }
            }
        }
        EXPECT_EQ(total_dists, 114);
    }

    // Sample 3: chain A, seqnum 4, PHE -> 209 total distances
    {
        auto it = res_lookup.find(res_key{"A", 4});
        ASSERT_NE(it, res_lookup.end());
        gemmi::CRA cra = it->second;
        EXPECT_EQ(cra.chain->name, "A");
        EXPECT_EQ(cra.residue->seqid.num.value, 4);
        EXPECT_EQ(cra.residue->name, "PHE");
        int total_dists = 0;
        for (auto const& outer : residue_distances_map) {
            if (outer.first.residue == cra.residue) {
                for (auto const& inner : outer.second) {
                    total_dists += static_cast<int>(inner.second.size());
                }
            }
        }
        EXPECT_EQ(total_dists, 209);
    }

    // --- Case: empty selection (should be a no-op) ---
    std::vector<gemmi::CRA> empty_selection;

    residue_distances_map.clear();
    before_generate_size = residue_distances_map.size();
    EXPECT_EQ(before_generate_size, 0u);

    coot::typed_distances::generate_gemmi(model, empty_selection, residue_distances_map, residues_within_window);

    after_generate_size = residue_distances_map.size();
    EXPECT_EQ(after_generate_size, 0u);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}