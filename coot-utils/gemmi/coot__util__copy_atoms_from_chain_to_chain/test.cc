#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, copy_atoms_from_chain_to_chain) {
    // Load the PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    // Get chain A from the structure
    gemmi::Model& model = st.models[0];
    gemmi::Chain* chain_A = nullptr;
    for (auto& chain : model.chains) {
        if (chain.name == "A") {
            chain_A = &chain;
            break;
        }
    }
    ASSERT_NE(chain_A, nullptr) << "Chain A not found";

    // Case 1: Valid copy between chains with matching structure
    {
        // Create a new chain and copy residues one by one
        gemmi::Chain chain_B_copy;
        chain_B_copy.name = "B";

        for (const auto& res : chain_A->residues) {
            chain_B_copy.residues.push_back(res);
        }

        EXPECT_EQ(chain_A->residues.size(), 267);
        EXPECT_EQ(chain_B_copy.residues.size(), 267);

        // Modify chain_B_copy coordinates to be different initially
        for (auto& res : chain_B_copy.residues) {
            for (auto& at : res.atoms) {
                at.pos.x += 10.0;
                at.pos.y += 10.0;
                at.pos.z += 10.0;
            }
        }

        // Call the function
        coot::util::copy_atoms_from_chain_to_chain_gemmi(chain_A, &chain_B_copy);

        // Verify coordinates were copied
        bool success = true;
        for (size_t ires = 0; ires < chain_A->residues.size(); ++ires) {
            const gemmi::Residue& res_A = chain_A->residues[ires];
            const gemmi::Residue& res_B = chain_B_copy.residues[ires];
            if (res_A.atoms.size() == res_B.atoms.size()) {
                for (size_t iat = 0; iat < res_A.atoms.size(); ++iat) {
                    const gemmi::Atom& at_A = res_A.atoms[iat];
                    const gemmi::Atom& at_B = res_B.atoms[iat];
                    double diff_sq = std::pow(at_A.pos.x - at_B.pos.x, 2) +
                                     std::pow(at_A.pos.y - at_B.pos.y, 2) +
                                     std::pow(at_A.pos.z - at_B.pos.z, 2);
                    if (std::sqrt(diff_sq) > 0.001) {
                        success = false;
                    }
                }
            }
        }
        EXPECT_TRUE(success);
    }

    // Case 2: Mismatched residue counts - should print error
    {
        gemmi::Chain short_chain;
        short_chain.name = "X";

        size_t n_add = std::min(size_t(10), chain_A->residues.size());
        for (size_t ires = 0; ires < n_add; ++ires) {
            short_chain.residues.push_back(chain_A->residues[ires]);
        }

        EXPECT_EQ(chain_A->residues.size(), 267);
        EXPECT_EQ(short_chain.residues.size(), 10);

        // This should print an error message
        coot::util::copy_atoms_from_chain_to_chain_gemmi(chain_A, &short_chain);
    }

    // Case 3: Mismatched atom counts within same residue - should print error
    {
        gemmi::Chain atom_short_chain;
        atom_short_chain.name = "Y";

        // Copy all residues
        for (const auto& res : chain_A->residues) {
            atom_short_chain.residues.push_back(res);
        }

        // Get first residue and delete some atoms
        if (!atom_short_chain.residues.empty()) {
            gemmi::Residue& res = atom_short_chain.residues[0];
            if (res.atoms.size() > 2) {
                // Delete atoms (keep only first 2)
                while (res.atoms.size() > 2) {
                    res.atoms.pop_back();
                }
            }
        }

        // This should print an error message about mismatched atom count
        coot::util::copy_atoms_from_chain_to_chain_gemmi(chain_A, &atom_short_chain);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}