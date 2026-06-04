#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, fill_ligand_atom_neighbour_map) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Model &model = st.models[0];

    // Get chain A, residue ILE at position 10
    gemmi::Residue *res = nullptr;
    for (gemmi::Chain &chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue &r : chain.residues) {
                if (r.name == "ILE" && r.seqid.num.value == 10) {
                    res = &r;
                    break;
                }
            }
        }
    }
    ASSERT_NE(res, nullptr);

    EXPECT_STREQ(res->name.c_str(), "ILE");
    EXPECT_EQ(res->seqid.num.value, 10);

    // Case 1: valid model — clear map from constructor init(), then call fill
    {
        coot::atom_overlaps_container_t aoc;
        aoc.model = &model;
        aoc.res_central = res;

        aoc.ligand_atom_neighbour_map.clear();

        int before = static_cast<int>(aoc.ligand_atom_neighbour_map.size());
        EXPECT_EQ(before, 0);

        aoc.fill_ligand_atom_neighbour_map_gemmi();

        int after = static_cast<int>(aoc.ligand_atom_neighbour_map.size());
        EXPECT_EQ(after, 8);

        // Check first 3 entries as oracle printed
        auto it = aoc.ligand_atom_neighbour_map.begin();
        EXPECT_NE(it, aoc.ligand_atom_neighbour_map.end());
        EXPECT_EQ(it->first, 0);
        EXPECT_EQ(static_cast<int>(it->second.size()), 1);

        ++it;
        EXPECT_NE(it, aoc.ligand_atom_neighbour_map.end());
        EXPECT_EQ(it->first, 1);
        EXPECT_EQ(static_cast<int>(it->second.size()), 3);

        ++it;
        EXPECT_NE(it, aoc.ligand_atom_neighbour_map.end());
        EXPECT_EQ(it->first, 2);
        EXPECT_EQ(static_cast<int>(it->second.size()), 2);
    }

    // Case 2: set model to nullptr on an already-constructed object (guard path)
    {
        coot::atom_overlaps_container_t aoc_null;
        aoc_null.model = &model;
        aoc_null.res_central = res;

        aoc_null.ligand_atom_neighbour_map.clear();
        aoc_null.model = nullptr; // trigger the if(mol) guard

        int before = static_cast<int>(aoc_null.ligand_atom_neighbour_map.size());
        EXPECT_EQ(before, 0);

        aoc_null.fill_ligand_atom_neighbour_map_gemmi();

        int after = static_cast<int>(aoc_null.ligand_atom_neighbour_map.size());
        EXPECT_EQ(after, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}