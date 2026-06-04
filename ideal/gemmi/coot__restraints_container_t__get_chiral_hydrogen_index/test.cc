#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_chiral_hydrogen_index) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    EXPECT_EQ(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];

    // Find chain A
    gemmi::Chain* chain = nullptr;
    for (gemmi::Chain& ch : model.chains) {
        if (ch.name == "A") {
            chain = &ch;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    // Collect all atoms from chain A in order, building is_hydrogen_flags
    std::vector<bool> is_hydrogen_flags;
    for (gemmi::Residue& res : chain->residues) {
        for (gemmi::Atom& atom : res.atoms) {
            is_hydrogen_flags.push_back(atom.element.is_hydrogen());
        }
    }
    ASSERT_FALSE(is_hydrogen_flags.empty());

    // Find a hydrogen and non-hydrogen atom index
    int hydrogen_index = -1;
    int non_hydrogen_index = -1;

    for (size_t i = 0; i < is_hydrogen_flags.size(); i++) {
        if (is_hydrogen_flags[i] && hydrogen_index < 0) {
            hydrogen_index = static_cast<int>(i);
        } else if (!is_hydrogen_flags[i] && non_hydrogen_index < 0) {
            non_hydrogen_index = static_cast<int>(i);
        }
        if (hydrogen_index >= 0 && non_hydrogen_index >= 0) break;
    }

    EXPECT_GE(hydrogen_index, 0);
    EXPECT_GE(non_hydrogen_index, 0);

    // Find a second non-hydrogen index (different from first)
    int non_hydrogen_index2 = -1;
    for (size_t i = 0; i < is_hydrogen_flags.size(); i++) {
        if (static_cast<int>(i) != non_hydrogen_index && !is_hydrogen_flags[i]) {
            non_hydrogen_index2 = static_cast<int>(i);
            break;
        }
    }
    EXPECT_GE(non_hydrogen_index2, 0);

    int indexc = non_hydrogen_index;

    // Case 1: exactly one hydrogen among index_1, index_2, index_3
    {
        int result = coot::get_chiral_hydrogen_index_gemmi(
            hydrogen_index, non_hydrogen_index, non_hydrogen_index2,
            is_hydrogen_flags);
        EXPECT_EQ(result, 8);
    }

    // Case 2: no hydrogens among the three → should return -1
    {
        int result = coot::get_chiral_hydrogen_index_gemmi(
            non_hydrogen_index, non_hydrogen_index2, non_hydrogen_index,
            is_hydrogen_flags);
        EXPECT_EQ(result, -1);
    }

    // Case 3: two hydrogens among the three → should return -1
    int hydrogen_index2 = -1;
    for (size_t i = 0; i < is_hydrogen_flags.size(); i++) {
        if (static_cast<int>(i) != hydrogen_index && is_hydrogen_flags[i]) {
            hydrogen_index2 = static_cast<int>(i);
            break;
        }
    }

    if (hydrogen_index2 >= 0) {
        int result = coot::get_chiral_hydrogen_index_gemmi(
            hydrogen_index, hydrogen_index2, non_hydrogen_index,
            is_hydrogen_flags);
        EXPECT_EQ(result, -1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}