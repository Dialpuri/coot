#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <utility>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "function.hh"

// Build typed_atoms from the model — a simplified version of the original's make_typed_atoms
// In the original, this uses protein_geometry to assign atom types.
// Here we use a simple approach: the "type" is the first letter of the element.
static std::vector<std::pair<gemmi::CRA, std::string>>
make_typed_atoms_gemmi(gemmi::Model& model) {
    std::vector<std::pair<gemmi::CRA, std::string>> result;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            for (gemmi::Atom& atom : res.atoms) {
                gemmi::CRA cra{&chain, &res, &atom};
                // Use the element name as the "type" (first char of element name)
                std::string type = atom.element.name();
                result.push_back({cra, type});
            }
        }
    }
    return result;
}

TEST(OracleTest, calculate_daca) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Gemmi models are 0-indexed
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // make typed atoms
    std::vector<std::pair<gemmi::CRA, std::string>> ta = make_typed_atoms_gemmi(model);

    // Find chain A
    gemmi::Chain* chain_a = nullptr;
    for (gemmi::Chain& ch : model.chains) {
        if (ch.name == std::string("A")) {
            chain_a = &ch;
            break;
        }
    }
    ASSERT_NE(chain_a, nullptr);

    // Case 1: REFERENCE mode on residue at index 10 (seqnum 11, GLY)
    {
        ASSERT_LT(10, chain_a->residues.size());
        gemmi::Residue& res = chain_a->residues[10];
        gemmi::CRA ref_cra{chain_a, &res, nullptr};

        int result = coot::daca::calculate_daca_gemmi(ref_cra, ta, coot::daca::REFERENCE);
        EXPECT_EQ(result, 0);
    }

    // Case 2: ANALYSIS mode on the same residue (after REFERENCE data built)
    {
        gemmi::Residue& res = chain_a->residues[10];
        gemmi::CRA ref_cra{chain_a, &res, nullptr};

        int result = coot::daca::calculate_daca_gemmi(ref_cra, ta, coot::daca::ANALYSIS);
        EXPECT_EQ(result, 0);
    }

    // Case 3: ANALYSIS mode on a different residue at index 50 (seqnum 60, HIS)
    {
        ASSERT_LT(50, chain_a->residues.size());
        gemmi::Residue& res = chain_a->residues[50];
        gemmi::CRA ref_cra{chain_a, &res, nullptr};

        int result = coot::daca::calculate_daca_gemmi(ref_cra, ta, coot::daca::ANALYSIS);
        EXPECT_EQ(result, 0);
    }

    // Case 4: REFERENCE mode on a different residue (HIS at index 50) — complementary case
    {
        gemmi::Residue& res = chain_a->residues[50];
        gemmi::CRA ref_cra{chain_a, &res, nullptr};

        int result = coot::daca::calculate_daca_gemmi(ref_cra, ta, coot::daca::REFERENCE);
        EXPECT_EQ(result, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}