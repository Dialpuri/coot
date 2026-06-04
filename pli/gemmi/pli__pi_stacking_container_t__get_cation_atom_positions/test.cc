#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include "function.hh"

// Helper to find a residue by chain name and sequence number in a gemmi Structure
static const gemmi::Residue* find_residue(const gemmi::Structure &st,
                                           char chain_id,
                                           int seq_num) {
    if (st.models.empty()) return nullptr;
    const gemmi::Model &model = st.models[0];
    for (const gemmi::Chain &chain : model.chains) {
        if (chain.name == std::string(1, chain_id)) {
            for (const gemmi::Residue &res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, get_cation_atom_positions) {
    gemmi::Structure st = gemmi::read_structure_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");

    pli::pi_stacking_container_t container;

    // Case 1: LYS residue (seq 6) - should find NZ atom
    {
        const gemmi::Residue *res_ptr = find_residue(st, 'A', 6);
        ASSERT_NE(res_ptr, nullptr);
        EXPECT_EQ(res_ptr->name, "LYS");
        std::vector<clipper::Coord_orth> positions =
            container.get_cation_atom_positions_gemmi(*res_ptr);
        EXPECT_EQ(positions.size(), 1u);
        EXPECT_NEAR(positions[0].x(), 17.011, 1e-3);
        EXPECT_NEAR(positions[0].y(), 11.836, 1e-3);
        EXPECT_NEAR(positions[0].z(), 73.255, 1e-3);
    }

    // Case 2: ARG residue (seq 22) - should find NH1 and NH2 atoms
    {
        const gemmi::Residue *res_ptr = find_residue(st, 'A', 22);
        ASSERT_NE(res_ptr, nullptr);
        EXPECT_EQ(res_ptr->name, "ARG");
        std::vector<clipper::Coord_orth> positions =
            container.get_cation_atom_positions_gemmi(*res_ptr);
        EXPECT_EQ(positions.size(), 2u);
        EXPECT_NEAR(positions[0].x(), 27.152, 1e-3);
        EXPECT_NEAR(positions[0].y(), 5.366, 1e-3);
        EXPECT_NEAR(positions[0].z(), 80.049, 1e-3);
        EXPECT_NEAR(positions[1].x(), 29.226, 1e-3);
        EXPECT_NEAR(positions[1].y(), 4.955, 1e-3);
        EXPECT_NEAR(positions[1].z(), 79.17, 1e-3);
    }

    // Case 3: Non-cation residue ALA (seq 21) - should return empty vector
    {
        const gemmi::Residue *res_ptr = find_residue(st, 'A', 21);
        ASSERT_NE(res_ptr, nullptr);
        EXPECT_EQ(res_ptr->name, "ALA");
        std::vector<clipper::Coord_orth> positions =
            container.get_cation_atom_positions_gemmi(*res_ptr);
        EXPECT_EQ(positions.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}