#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <map>
#include <set>
#include <vector>

TEST(OracleTest, add_fixed_atoms_from_non_bonded_neighbours) {
    // Load the test structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model &model = st.models[0];

    // Build atom index map: pointer -> sequential global index
    std::map<gemmi::Atom*, int> atom_index_map;
    int idx = 0;
    for (gemmi::Chain &chain : model.chains) {
        for (gemmi::Residue &residue : chain.residues) {
            for (gemmi::Atom &atom : residue.atoms) {
                atom_index_map[&atom] = idx++;
            }
        }
    }

    // Find non-bonded neighbour residues: A/98 GLY, A/102 PRO, A/103 LEU, A/104 ILE, A/105 LYS
    std::vector<gemmi::Residue*> non_bonded_neighbour_residues;
    for (gemmi::Chain &chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue &residue : chain.residues) {
                int seqnum = residue.seqid.num.value;
                if (seqnum == 98 || seqnum == 102 || seqnum == 103 || seqnum == 104 || seqnum == 105) {
                    non_bonded_neighbour_residues.push_back(&residue);
                }
            }
        }
    }

    EXPECT_EQ(non_bonded_neighbour_residues.size(), std::size_t(5));

    // Case 1: Clear fixed_atom_indices, then re-call to see it repopulate
    std::set<int> fixed_atom_indices;
    std::size_t fixed_before = fixed_atom_indices.size();
    EXPECT_EQ(fixed_before, std::size_t(0));

    add_fixed_atoms_from_non_bonded_neighbours_gemmi(non_bonded_neighbour_residues, fixed_atom_indices, atom_index_map);

    std::size_t fixed_after = fixed_atom_indices.size();
    EXPECT_EQ(fixed_after, std::size_t(36));

    // Check non-bonded neighbour residues
    EXPECT_EQ(non_bonded_neighbour_residues.size(), std::size_t(5));

    {
        gemmi::Residue *r = non_bonded_neighbour_residues[0];
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->seqid.num.value, 98);
        EXPECT_EQ(r->name, "GLY");
    }
    {
        gemmi::Residue *r = non_bonded_neighbour_residues[1];
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->seqid.num.value, 102);
        EXPECT_EQ(r->name, "PRO");
    }
    {
        gemmi::Residue *r = non_bonded_neighbour_residues[2];
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->seqid.num.value, 103);
        EXPECT_EQ(r->name, "LEU");
    }
    {
        gemmi::Residue *r = non_bonded_neighbour_residues[3];
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->seqid.num.value, 104);
        EXPECT_EQ(r->name, "ILE");
    }
    {
        gemmi::Residue *r = non_bonded_neighbour_residues[4];
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->seqid.num.value, 105);
        EXPECT_EQ(r->name, "LYS");
    }

    // Case 2: Call again with non-empty fixed_atom_indices (idempotent)
    std::size_t fixed_before_2 = fixed_atom_indices.size();
    EXPECT_EQ(fixed_before_2, std::size_t(36));

    add_fixed_atoms_from_non_bonded_neighbours_gemmi(non_bonded_neighbour_residues, fixed_atom_indices, atom_index_map);

    std::size_t fixed_after_2 = fixed_atom_indices.size();
    EXPECT_EQ(fixed_after_2, std::size_t(36));
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}