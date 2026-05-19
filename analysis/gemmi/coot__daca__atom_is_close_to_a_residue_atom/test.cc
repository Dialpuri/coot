#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, daca_atom_is_close_to_a_residue_atom) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get a reference residue (e.g., LYS 10 in chain A)
    gemmi::Residue* ref_res = nullptr;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (std::string(1, res.seqid.icode) == " " && res.seqid.num.value == 10) {
                    ref_res = &res;
                    break;
                }
            }
            if (ref_res) break;
        }
        if (ref_res) break;
    }
    
    ASSERT_NE(ref_res, nullptr) << "Failed to get residue 10 in chain A";
    
    // Case 1: atom close to the reference residue
    gemmi::Atom test_atom_close;
    if (!ref_res->atoms.empty()) {
        const gemmi::Atom& ref_atom = ref_res->atoms[0];
        test_atom_close.pos.x = ref_atom.pos.x + 0.5;
        test_atom_close.pos.y = ref_atom.pos.y + 0.5;
        test_atom_close.pos.z = ref_atom.pos.z + 0.5;
        test_atom_close.name = "CA";
    } else {
        // Fallback: use a dummy position near the residue's center
        test_atom_close.pos.x = ref_res->atoms.front().pos.x + 0.5;
        test_atom_close.pos.y = ref_res->atoms.front().pos.y + 0.5;
        test_atom_close.pos.z = ref_res->atoms.front().pos.z + 0.5;
        test_atom_close.name = "CA";
    }
    
    bool is_close = coot::daca::atom_is_close_to_a_residue_atom_gemmi(test_atom_close, *ref_res);
    EXPECT_TRUE(is_close);
    
    // Case 2: atom far from the reference residue
    gemmi::Atom test_atom_far;
    test_atom_far.pos.x = 100.0;
    test_atom_far.pos.y = 100.0;
    test_atom_far.pos.z = 100.0;
    test_atom_far.name = "CA";
    
    bool is_far = coot::daca::atom_is_close_to_a_residue_atom_gemmi(test_atom_far, *ref_res);
    EXPECT_FALSE(is_far);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}