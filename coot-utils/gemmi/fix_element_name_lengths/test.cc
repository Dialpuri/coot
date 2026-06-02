#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, fix_element_name_lengths) {
    // Case from oracle: basic case with example.pdb
    {
        gemmi::Structure st = gemmi::read_pdb_file(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        auto& model = st.models[0];
        auto& chain = model.chains[0];
        auto& res = chain.residues[0];
        auto& first_atom = res.atoms[0];

        std::string elem_before = first_atom.element.name();

        int single_count_before = 0;
        for (auto& m : st.models)
            for (auto& ch : m.chains)
                for (auto& r : ch.residues)
                    for (auto& a : r.atoms) {
                        std::string en = a.element.name();
                        if (en.length() == 1)
                            single_count_before++;
                    }

        EXPECT_EQ(first_atom.name, "N");
        EXPECT_EQ(elem_before, "N");
        // In gemmi element names are unpadded, so single_count reflects real
        // single-char elements (C, N, O, H, S, etc.), not a padding artifact.

        fix_element_name_lengths_gemmi(st);

        std::string elem_after = first_atom.element.name();

        int single_count_after = 0;
        for (auto& m : st.models)
            for (auto& ch : m.chains)
                for (auto& r : ch.residues)
                    for (auto& a : r.atoms) {
                        std::string en = a.element.name();
                        if (en.length() == 1)
                            single_count_after++;
                    }

        EXPECT_EQ(elem_after, "N");
        EXPECT_EQ(single_count_after, single_count_before);
    }

    // Case 2: element remains unchanged after fix (gemmi has no padding)
    {
        gemmi::Structure st = gemmi::read_pdb_file(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        auto& model = st.models[0];
        auto& chain = model.chains[0];
        auto& res = chain.residues[0];
        auto& test_atom = res.atoms[0];

        test_atom.element = gemmi::El::N;

        std::string ele_before = test_atom.element.name();
        EXPECT_EQ(ele_before.length(), 1u);

        fix_element_name_lengths_gemmi(st);

        std::string ele_after = test_atom.element.name();
        EXPECT_EQ(ele_after.length(), 1u);
        EXPECT_EQ(ele_after, "N");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}