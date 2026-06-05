#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, restraints_container_t_add_link_bond_gemmi) {
    // Load the PDB and get the model
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    EXPECT_GT(st.models.size(), 0u);
    gemmi::Model& model = st.models[0];

    coot::protein_geometry pg;
    pg.init_standard();

    // Get two adjacent residues from chain A: residue 1 and residue 2
    const gemmi::Chain* chain = nullptr;
    for (gemmi::Chain& ch : model.chains) {
        if (ch.name == "A") { chain = &ch; break; }
    }
    ASSERT_NE(chain, nullptr);

    // Iterate over residues to find ones with seqnum 1 and 2
    gemmi::Residue const *res1 = nullptr;
    gemmi::Residue const *res2 = nullptr;

    for (const gemmi::Residue& res : chain->residues) {
        if (res.seqid.num.value == 1 && !res1) res1 = &res;
        if (res.seqid.num.value == 2 && !res2) res2 = &res;
    }

    EXPECT_NE(res1, nullptr);
    EXPECT_NE(res2, nullptr);

    // Prepare the data structures used by add_link_bond_gemmi
    std::vector<coot::simple_restraint> restraints_vec;
    std::vector<std::set<int>> bonded_atom_indices;

    // Pre-size bonded_atom_indices to hold at least 100 atoms per residue
    // (index = residue_index * 100 + atom_index)
    bonded_atom_indices.resize(300);

    // Case 1: Valid TRANS link between adjacent residues 1 and 2
    {
        std::string link_type = "TRANS";
        short int is_fixed_first = 0;
        short int is_fixed_second = 0;

        int bonds_before = static_cast<int>(restraints_vec.size());
        int nbond = coot::add_link_bond_gemmi(
            restraints_vec, bonded_atom_indices,
            link_type, res1, res2,
            is_fixed_first, is_fixed_second,
            pg, 0, 1); // residue indices 0 and 1
        int bonds_after = static_cast<int>(restraints_vec.size());

        EXPECT_EQ(nbond, 1);
        EXPECT_EQ(bonds_before, 0);
        EXPECT_EQ(bonds_after, 1);
    }

    // Case 2: Invalid link type — should produce 0 bonds
    {
        std::string link_type = "NONEXISTENT";
        short int is_fixed_first = 0;
        short int is_fixed_second = 0;

        int bonds_before = static_cast<int>(restraints_vec.size());
        int nbond = coot::add_link_bond_gemmi(
            restraints_vec, bonded_atom_indices,
            link_type, res1, res2,
            is_fixed_first, is_fixed_second,
            pg, 0, 1);
        int bonds_after = static_cast<int>(restraints_vec.size());

        EXPECT_EQ(nbond, 0);
        EXPECT_EQ(bonds_before, 1);
        EXPECT_EQ(bonds_after, 1);
    }

    // Case 3: Additional complementary case — valid TRANS link between
    // residues 2 and 3 (different residues, to test mutation on a
    // non-empty restraints_vec)
    {
        gemmi::Residue const *res3 = nullptr;
        for (const gemmi::Residue& res : chain->residues) {
            if (res.seqid.num.value == 3 && !res3) res3 = &res;
        }
        if (res3) {
            std::string link_type = "TRANS";
            short int is_fixed_first = 0;
            short int is_fixed_second = 0;

            int bonds_before = static_cast<int>(restraints_vec.size());
            int nbond = coot::add_link_bond_gemmi(
                restraints_vec, bonded_atom_indices,
                link_type, res2, res3,
                is_fixed_first, is_fixed_second,
                pg, 1, 2);
            int bonds_after = static_cast<int>(restraints_vec.size());

            EXPECT_GT(nbond, 0);
            EXPECT_EQ(bonds_after, bonds_before + nbond);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}