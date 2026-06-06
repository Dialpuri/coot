#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, calculate_residue_normals) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    if (st.models.empty()) {
        FAIL() << "No models loaded";
    }

    gemmi::Model& model = st.models[0];

    // Chain A (first chain)
    if (model.chains.empty()) {
        FAIL() << "No chains found";
    }

    gemmi::Chain& chain = model.chains[0];

    // Case 1: Atoms from first residue (MET 1) — protein residue, no base normal expected
    {
        gemmi::Residue& res = chain.residues[0];
        EXPECT_EQ(res.name, "MET");

        std::vector<const gemmi::Residue*> residues_vec;
        residues_vec.push_back(&res);

        std::map<const gemmi::Residue*, clipper::Coord_orth> normals =
            coot::stack_and_pair::calculate_residue_normals_gemmi(residues_vec);

        EXPECT_EQ(normals.size(), 0);
    }

    // Case 2: Atoms from a different protein residue (PHE 4) — also no base normal
    {
        gemmi::Residue& res = chain.residues[3];
        EXPECT_EQ(res.name, "PHE");

        std::vector<const gemmi::Residue*> residues_vec;
        residues_vec.push_back(&res);

        std::map<const gemmi::Residue*, clipper::Coord_orth> normals =
            coot::stack_and_pair::calculate_residue_normals_gemmi(residues_vec);

        EXPECT_EQ(normals.size(), 0);
    }

    // Case 3: Empty atom list (edge case)
    {
        std::vector<const gemmi::Residue*> residues_vec;

        std::map<const gemmi::Residue*, clipper::Coord_orth> normals =
            coot::stack_and_pair::calculate_residue_normals_gemmi(residues_vec);

        EXPECT_EQ(normals.size(), 0);
    }

    // Case 4: Multiple residues — tests dedup via done_res set
    {
        std::vector<const gemmi::Residue*> residues_vec;
        // Grab first three residues (or fewer if chain is short)
        int n_res = std::min<int>(chain.residues.size(), 3);
        for (int r = 0; r < n_res; r++) {
            residues_vec.push_back(&chain.residues[r]);
        }

        std::map<const gemmi::Residue*, clipper::Coord_orth> normals =
            coot::stack_and_pair::calculate_residue_normals_gemmi(residues_vec);

        EXPECT_EQ(normals.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}