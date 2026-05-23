#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper function to pad atom names as MMDB does
std::string mmdb_atom_name(const std::string& name) {
    if (name.length() == 4) {
        return name;
    } else if (name.length() < 4) {
        std::string result = " " + name;
        while (result.length() < 4) {
            result += " ";
        }
        return result;
    } else {
        return name.substr(0, 4);
    }
}

TEST(OracleTest, get_quads_ALA_residue_21) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get residue A/21 (ALA) which has atoms: C, CA, CB, N, O
    gemmi::Model& model = st.models[0];
    gemmi::Residue* res = nullptr;
    for (auto& chain : model.chains) {
        for (auto& r : chain.residues) {
            if (r.name == "ALA" && r.seqid.num.value == 21) {
                res = &r;
                break;
            }
        }
        if (res) break;
    }
    ASSERT_NE(res, nullptr) << "Could not get residue A/21";

    std::vector<coot::monomer_utils::atom_name_quad> quads;

    // Case 1: Valid quad - ALA has C, CA, CB, N
    quads.push_back(coot::monomer_utils::atom_name_quad(mmdb_atom_name("C"), mmdb_atom_name("CA"), mmdb_atom_name("CB"), mmdb_atom_name("N")));
    
    // Case 2: Another valid quad - ALA has CA, CB, N, O
    quads.push_back(coot::monomer_utils::atom_name_quad(mmdb_atom_name("CA"), mmdb_atom_name("CB"), mmdb_atom_name("N"), mmdb_atom_name("O")));
    
    // Case 3: Invalid quad with non-existent atoms
    quads.push_back(coot::monomer_utils::atom_name_quad(mmdb_atom_name("XX"), mmdb_atom_name("YY"), mmdb_atom_name("ZZ"), mmdb_atom_name("AA")));

    // Call get_quads_gemmi
    std::vector<coot::monomer_utils::atom_index_quad> result = coot::monomer_utils::get_quads_gemmi(quads, *res);

    // Expected: only 2 valid quads found (the 3rd has non-existent atoms)
    EXPECT_EQ(result.size(), 2);

    // First quad: C, CA, CB, N -> indices (2, 1, 4, 0) based on output
    EXPECT_EQ(result[0].index1, 2);
    EXPECT_EQ(result[0].index2, 1);
    EXPECT_EQ(result[0].index3, 4);
    EXPECT_EQ(result[0].index4, 0);

    // Second quad: CA, CB, N, O -> indices (1, 4, 0, 3) based on output
    EXPECT_EQ(result[1].index1, 1);
    EXPECT_EQ(result[1].index2, 4);
    EXPECT_EQ(result[1].index3, 0);
    EXPECT_EQ(result[1].index4, 3);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}