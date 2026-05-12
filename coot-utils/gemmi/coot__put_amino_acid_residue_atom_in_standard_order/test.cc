#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper function to get atom names from a residue as a string
std::string get_atom_names(gemmi::Residue& res) {
    std::string result;
    for (size_t i = 0; i < res.atoms.size(); i++) {
        result += res.atoms[i].name;
        if (i < res.atoms.size() - 1) result += ", ";
    }
    return result;
}

// Helper function to get the first atom name from a residue
std::string get_first_atom_name(gemmi::Residue& res) {
    if (res.atoms.size() > 0) {
        return res.atoms[0].name;
    }
    return "";
}

TEST(OracleTest, PutAminoAcidResidueAtomInStandardOrder) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB file";

    // Test 1: Normal residue (GLU 2) - N is already first
    {
        gemmi::Model& model = st.models[0];
        bool found = false;
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 2) {
                    found = true;
                    std::string input_atoms = get_atom_names(res);
                    coot::put_amino_acid_residue_atom_in_standard_order_gemmi(res);
                    std::string output_atoms = get_atom_names(res);
                    std::string first_atom = get_first_atom_name(res);
                    
                    EXPECT_EQ(input_atoms, "N, CA, C, O, CB, CG, CD, OE1, OE2");
                    EXPECT_EQ(output_atoms, "N, CA, C, O, CB, CG, CD, OE1, OE2");
                    EXPECT_EQ(first_atom, "N");
                    break;
                }
            }
            if (found) break;
        }
        ASSERT_TRUE(found) << "Residue 2 not found";
    }

    // Test 2: GLN 5
    {
        gemmi::Model& model = st.models[0];
        bool found = false;
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 5) {
                    found = true;
                    std::string input_atoms = get_atom_names(res);
                    coot::put_amino_acid_residue_atom_in_standard_order_gemmi(res);
                    std::string output_atoms = get_atom_names(res);
                    std::string first_atom = get_first_atom_name(res);
                    
                    EXPECT_EQ(input_atoms, "N, CA, C, O, CB, CG, CD, OE1, NE2");
                    EXPECT_EQ(output_atoms, "N, CA, C, O, CB, CG, CD, OE1, NE2");
                    EXPECT_EQ(first_atom, "N");
                    break;
                }
            }
            if (found) break;
        }
        ASSERT_TRUE(found) << "Residue 5 not found";
    }

    // Test 3: LYS 33
    {
        gemmi::Model& model = st.models[0];
        bool found = false;
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 33) {
                    found = true;
                    std::string input_atoms = get_atom_names(res);
                    coot::put_amino_acid_residue_atom_in_standard_order_gemmi(res);
                    std::string output_atoms = get_atom_names(res);
                    std::string first_atom = get_first_atom_name(res);
                    
                    EXPECT_EQ(input_atoms, "N, CA, C, O, CB, CG, CD, CE, NZ");
                    EXPECT_EQ(output_atoms, "N, CA, C, O, CB, CG, CD, CE, NZ");
                    EXPECT_EQ(first_atom, "N");
                    break;
                }
            }
            if (found) break;
        }
        ASSERT_TRUE(found) << "Residue 33 not found";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}