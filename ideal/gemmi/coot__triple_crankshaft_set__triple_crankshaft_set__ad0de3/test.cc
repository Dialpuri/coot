#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, TripleCrankshaftSetConstructor) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get 6 consecutive residues from chain A
    gemmi::Residue *res[6];
    const char* seq_ids[] = {"10", "11", "12", "13", "14", "15"};
    
    for (int i = 0; i < 6; i++) {
        bool found = false;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& residue : chain.residues) {
                    if (std::to_string(residue.seqid.num.value) == seq_ids[i]) {
                        res[i] = &residue;
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (found) break;
        }
        ASSERT_NE(res[i], nullptr) << "Failed to get residue: " << seq_ids[i];
    }

    // Build residue types vector
    std::vector<std::string> residue_types;
    for (int i = 0; i < 6; i++) {
        residue_types.push_back(res[i]->name);
    }

    // Test case 1: valid constructor call with consecutive residues
    {
        // Print input for debugging
        std::cout << "INPUT  res[0]: " << seq_ids[0] << std::endl;
        std::cout << "INPUT  res[1]: " << seq_ids[1] << std::endl;
        std::cout << "INPUT  res[2]: " << seq_ids[2] << std::endl;
        std::cout << "INPUT  res[3]: " << seq_ids[3] << std::endl;
        std::cout << "INPUT  res[4]: " << seq_ids[4] << std::endl;
        std::cout << "INPUT  res[5]: " << seq_ids[5] << std::endl;
        std::cout << "INPUT  residue_types: ";
        for (size_t i = 0; i < residue_types.size(); i++) {
            if (i > 0) std::cout << ",";
            std::cout << residue_types[i];
        }
        std::cout << std::endl;

        // This should not throw
        EXPECT_NO_THROW(coot::triple_crankshaft_set tcs(res[0], res[1], res[2], res[3], res[4], res[5], residue_types));

        std::cout << "OUTPUT constructor_success: true" << std::endl;
    }

    // Test case 2: non-consecutive residues (edge case)
    const char* edge_seq_ids[] = {"10", "12", "14", "16", "18", "20"};
    gemmi::Residue *edge_res[6];
    bool all_found = true;
    for (int i = 0; i < 6; i++) {
        bool found = false;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& residue : chain.residues) {
                    if (std::to_string(residue.seqid.num.value) == edge_seq_ids[i]) {
                        edge_res[i] = &residue;
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (found) break;
        }
        if (!found) {
            all_found = false;
            break;
        }
    }

    if (all_found) {
        std::vector<std::string> edge_types;
        for (int i = 0; i < 6; i++) {
            edge_types.push_back(edge_res[i]->name);
        }

        std::cout << "INPUT  edge_res[0]: " << edge_seq_ids[0] << std::endl;
        std::cout << "INPUT  edge_res[1]: " << edge_seq_ids[1] << std::endl;
        std::cout << "INPUT  edge_res[2]: " << edge_seq_ids[2] << std::endl;
        std::cout << "INPUT  edge_res[3]: " << edge_seq_ids[3] << std::endl;
        std::cout << "INPUT  edge_res[4]: " << edge_seq_ids[4] << std::endl;
        std::cout << "INPUT  edge_res[5]: " << edge_seq_ids[5] << std::endl;

        // This should also not throw
        EXPECT_NO_THROW(coot::triple_crankshaft_set edge_tcs(edge_res[0], edge_res[1], edge_res[2], edge_res[3], edge_res[4], edge_res[5], edge_types));
        std::cout << "OUTPUT edge_constructor_success: true" << std::endl;
    } else {
        std::cout << "INPUT  edge_case_skipped: true" << std::endl;
        std::cout << "OUTPUT edge_case_reason: not all residues found" << std::endl;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}