#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atom_tree_t) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_GE(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];

    // Find chain "A", residue 10 (ILE)
    gemmi::Residue* find_residue = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 10) {
                    find_residue = &res;
                    break;
                }
            }
        }
    }
    EXPECT_NE(find_residue, nullptr);

    // Case 1: Valid residue with sensible contact indices
    {
        EXPECT_NE(find_residue, nullptr);

        std::string res_name = find_residue->name;
        EXPECT_EQ(res_name, "ILE");
        EXPECT_EQ(find_residue->seqid.num.value, 10);

        std::string altconf = "";

        // Build a linear contact graph: 0-1-2-3-4 (5 atoms)
        std::vector<std::vector<int>> contact_indices(5);
        contact_indices[0].push_back(1);
        contact_indices[1].push_back(0); contact_indices[1].push_back(2);
        contact_indices[2].push_back(1); contact_indices[2].push_back(3);
        contact_indices[3].push_back(2); contact_indices[3].push_back(4);
        contact_indices[4].push_back(3);

        int base_atom_index = 0;

        std::map<std::string, int> name_to_index;
        std::vector<coot::atom_vertex> atom_vertex_vec;

        coot::atom_tree_t_gemmi(contact_indices, base_atom_index, find_residue, altconf, name_to_index, atom_vertex_vec);
        EXPECT_EQ(name_to_index.size(), 8u);
        EXPECT_EQ(atom_vertex_vec.size(), 8u);
    }

    // Case 2: Null residue - should throw
    {
        gemmi::Residue* null_res = nullptr;
        std::vector<std::vector<int>> contact_indices;
        int base_atom_index = 0;
        std::string altconf = "";
        std::map<std::string, int> name_to_index;
        std::vector<coot::atom_vertex> atom_vertex_vec;

        EXPECT_THROW(
            coot::atom_tree_t_gemmi(contact_indices, base_atom_index, null_res, altconf, name_to_index, atom_vertex_vec),
            std::runtime_error
        );
    }

    // Case 3: Valid residue with empty contact_indices (single atom, no connections)
    {
        EXPECT_NE(find_residue, nullptr);
        std::string altconf = "";

        std::vector<std::vector<int>> contact_indices(1); // single atom, no connections
        int base_atom_index = 0;

        std::map<std::string, int> name_to_index;
        std::vector<coot::atom_vertex> atom_vertex_vec;

        coot::atom_tree_t_gemmi(contact_indices, base_atom_index, find_residue, altconf, name_to_index, atom_vertex_vec);
        EXPECT_EQ(name_to_index.size(), 8u);
        EXPECT_EQ(atom_vertex_vec.size(), 8u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}