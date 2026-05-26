#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, beam_in_linked_residue_get_residue) {
    // case: setup - load test PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case: get reference residue (ILE at chain A, residue 10)
    gemmi::Residue* res = nullptr;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& r : chain.residues) {
                if (r.seqid.num.value == 10 && chain.name == "A") {
                    res = &r;
                    break;
                }
            }
            if (res) break;
        }
        if (res) break;
    }
    ASSERT_NE(res, nullptr) << "Reference residue not found";
    
    std::string ref_res_name = res->name;
    EXPECT_EQ(ref_res_name, "ILE") << "Expected ILE residue";
    
    // case: test get_residue_gemmi with matching comp_id
    std::string comp_id = ref_res_name;
    gemmi::Residue* found_res = coot::beam_in_linked_residue::get_residue_gemmi(comp_id, st);
    
    EXPECT_NE(found_res, nullptr) << "Residue should be found for matching comp_id";
    
    if (found_res) {
        EXPECT_EQ(found_res->name, "ILE") << "Residue name should match";
        EXPECT_EQ(found_res->seqid.num.value, 10) << "Residue sequence number should be 10";
        EXPECT_EQ(found_res, res) << "Should return same residue object";
    }
    
    // case: test get_residue_gemmi with non-existent comp_id
    std::string non_existent = "XXX";
    gemmi::Residue* not_found = coot::beam_in_linked_residue::get_residue_gemmi(non_existent, st);
    
    EXPECT_EQ(not_found, nullptr) << "Residue should not be found for non-existent comp_id";
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}