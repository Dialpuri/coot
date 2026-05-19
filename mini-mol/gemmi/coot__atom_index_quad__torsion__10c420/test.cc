#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, AtomIndexQuadTorsion) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get residue ILE at position 10 in chain A
    const gemmi::Model& model = st.models[0];
    const gemmi::Chain* chain_a = nullptr;
    for (const auto& chain : model.chains) {
        if (chain.name == "A") {
            chain_a = &chain;
            break;
        }
    }
    ASSERT_NE(chain_a, nullptr) << "Chain A not found";
    
    const gemmi::Residue* res = nullptr;
    for (const auto& r : chain_a->residues) {
        if (r.seqid.num.value == 10) {
            res = &r;
            break;
        }
    }
    ASSERT_NE(res, nullptr) << "Failed to get residue //A/10";

    // Case 1: Test with valid indices (first 4 atoms)
    {
        EXPECT_EQ(res->name, "ILE");
        EXPECT_EQ(res->atoms.size(), 8);

        if (res->atoms.size() >= 4) {
            coot::atom_index_quad quad(0, 1, 2, 3);
            // Use EXPECT_NO_THROW for the valid case, then check the value
            EXPECT_NO_THROW(quad.torsion_gemmi(*res));
            double angle = quad.torsion_gemmi(*res);
            // The oracle printed: 174.318 (rounded)
            // Actual value is 174.31781459093162
            EXPECT_NEAR(angle, 174.318, 2e-4);
        } else {
            // Fallback if residue has fewer than 4 atoms
            // (but we know ILE has 8 atoms)
            FAIL() << "Residue has fewer than 4 atoms";
        }
    }

    // Case 2: Test with invalid indices (out of range)
    {
        coot::atom_index_quad quad(999, 998, 997, 996);
        // The oracle printed exception: "bad atom indexing in atom_index_quad::torsion()"
        EXPECT_THROW(quad.torsion_gemmi(*res), std::runtime_error);
        try {
            quad.torsion_gemmi(*res);
            FAIL() << "Expected std::runtime_error for invalid indices";
        } catch (const std::runtime_error &e) {
            EXPECT_NE(std::string(e.what()).find("bad atom indexing in atom_index_quad::torsion()"), 
                      std::string::npos);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}