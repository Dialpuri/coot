#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, TorsionGeneralAtomIndex) {
    // Load PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB file";

    // Get residue //A/10 (chain A, residue number 10, no insertion code)
    gemmi::Residue* res = nullptr;
    for (auto& model : st.models)
        for (auto& chain : model.chains)
            for (auto& r : chain.residues)
                if (chain.name == "A" && r.seqid.num.value == 10 && r.seqid.icode == ' ') {
                    res = &r;
                    break;
                }
    ASSERT_NE(res, nullptr) << "Failed to get residue //A/10";

    // Debug: print atom names
    for (size_t i = 0; i < res->atoms.size(); ++i) {
        fprintf(stderr, "Atom %zu: name='%s'\n", i, res->atoms[i].name.c_str());
    }

    // Test case 1: valid atom spec that should match (CA)
    {
        coot::atom_spec_t spec("A", 10, "", " CA ", "");
        int idx = coot::torsion_general::atom_index_gemmi(spec, *res);
        EXPECT_EQ(idx, 1) << "Expected CA atom index to be 1";
        if (idx >= 0) {
            std::string atom_name = res->atoms[idx].name;
            EXPECT_EQ(atom_name, "CA") << "Expected atom name to be 'CA'";  // gemmi stores without spaces
        }
    }

    // Test case 2: valid atom spec for different atom (CB)
    {
        coot::atom_spec_t spec("A", 10, "", " CB ", "");
        int idx = coot::torsion_general::atom_index_gemmi(spec, *res);
        EXPECT_EQ(idx, 4) << "Expected CB atom index to be 4";
    }

    // Test case 3: invalid atom spec that should not match (ZZZ)
    {
        coot::atom_spec_t spec("A", 10, "", "ZZZ", "");
        int idx = coot::torsion_general::atom_index_gemmi(spec, *res);
        EXPECT_EQ(idx, -1) << "Expected ZZZ atom index to be -1 (not found)";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}