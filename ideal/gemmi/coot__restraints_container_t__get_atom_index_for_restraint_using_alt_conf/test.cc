#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_atom_index_for_restraint_using_alt_conf) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Model &model = st.models[0];

    // Find chain "A"
    const gemmi::Chain *chain = model.find_chain("A");
    ASSERT_NE(chain, nullptr);

    // Find residue with seqid 10
    const gemmi::Residue *res = nullptr;
    for (const auto &r : chain->residues) {
        if (r.seqid.num.value == 10) {
            res = &r;
            break;
        }
    }
    ASSERT_NE(res, nullptr);

    // Case 1: look up CA atom with empty alt_conf
    {
        std::string atom_name = "CA";
        std::string alt_conf = "";
        int idx = coot::get_atom_index_for_restraint_using_alt_conf_gemmi(
            atom_name, alt_conf, *res);
        EXPECT_EQ(idx, -1);
    }

    // Case 2: look up non-existent atom name
    {
        std::string atom_name = "ZZZ";
        std::string alt_conf = "";
        int idx = coot::get_atom_index_for_restraint_using_alt_conf_gemmi(
            atom_name, alt_conf, *res);
        EXPECT_EQ(idx, -1);
    }

    // Case 3: look up N atom with wrong alt_conf
    {
        std::string atom_name = "N";
        std::string alt_conf = "A";
        int idx = coot::get_atom_index_for_restraint_using_alt_conf_gemmi(
            atom_name, alt_conf, *res);
        EXPECT_EQ(idx, -1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}