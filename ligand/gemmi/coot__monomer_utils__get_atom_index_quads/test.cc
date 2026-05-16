#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include "function.hh"

TEST(OracleTest, get_atom_index_quads_gemmi) {
    // Load the PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_FALSE(st.models.empty());

    // Get the first model
    gemmi::Model& model = st.models[0];

    // Traverse chains to find ALA 21
    gemmi::Residue* target_res = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            if (res.name == "ALA" && res.seqid.num.value == 21) {
                target_res = &res;
                break;
            }
        }
        if (target_res) break;
    }

    EXPECT_NE(target_res, nullptr);

    // Get atom names from the residue
    // MMDB stores atom names padded with leading spaces (e.g., " C", " CA", " N")
    // gemmi trims them, so we pad to match original MMDB behavior
    std::vector<std::string> atom_names;
    for (gemmi::Atom& atom : target_res->atoms) {
        std::string padded(4 - atom.name.size(), ' ');
        padded += atom.name;
        atom_names.push_back(padded);
    }

    // Case 1: valid quad - ALA has C, CA, CB, N (trimmed names)
    {
        std::vector<coot::atom_name_quad> quads;
        quads.emplace_back("C", "CA", "CB", "N");

        std::vector<coot::atom_index_quad> result =
            coot::monomer_utils::get_atom_index_quads_gemmi(quads, atom_names);

        // All quads fail to match (oracle ground truth: quads_found: 0)
        EXPECT_EQ(result.size(), 0);
    }

    // Case 2: invalid quad with nonexistent atom name
    {
        std::vector<coot::atom_name_quad> quads;
        quads.emplace_back("ZZZZ", "CA", "CB", "N");

        std::vector<coot::atom_index_quad> result =
            coot::monomer_utils::get_atom_index_quads_gemmi(quads, atom_names);

        EXPECT_EQ(result.size(), 0);
    }

    // Case 3: multiple quads at once
    {
        std::vector<coot::atom_name_quad> quads;
        quads.emplace_back("C", "CA", "CB", "N");
        quads.emplace_back("CA", "CB", "N", "O");

        std::vector<coot::atom_index_quad> result =
            coot::monomer_utils::get_atom_index_quads_gemmi(quads, atom_names);

        EXPECT_EQ(result.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}