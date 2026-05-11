#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, filter_residues_by_solvent_contact) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "PDB load failed";

    // Get a reference residue - first try 100, then 10
    const gemmi::Residue* res_ref = nullptr;
    const gemmi::Chain* ref_chain = nullptr;
    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 100) {
                    res_ref = &res;
                    ref_chain = &chain;
                    break;
                }
            }
            if (res_ref) break;
        }
        if (res_ref) break;
    }
    if (!res_ref) {
        for (const gemmi::Model& model : st.models) {
            for (const gemmi::Chain& chain : model.chains) {
                for (const gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 10) {
                        res_ref = &res;
                        ref_chain = &chain;
                        break;
                    }
                }
                if (res_ref) break;
            }
            if (res_ref) break;
        }
    }
    ASSERT_NE(res_ref, nullptr) << "No residue found to use as reference";

    // Get residues near the reference (manually scan all residues within 5.0A)
    // Exclude the reference residue itself from the list
    std::vector<const gemmi::Residue*> residues;
    const double search_radius = 5.0;
    const double search_radius_sq = search_radius * search_radius;

    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& res : chain.residues) {
                // Skip the reference residue itself
                if (&res == res_ref) continue;
                
                // Check if this residue is within search_radius of res_ref
                bool within_range = false;
                for (const gemmi::Atom& ref_atom : res_ref->atoms) {
                    if (within_range) break;
                    for (const gemmi::Atom& atom : res.atoms) {
                        double dist_sq = (ref_atom.pos - atom.pos).length_sq();
                        if (dist_sq < search_radius_sq) {
                            within_range = true;
                            break;
                        }
                    }
                }
                if (within_range) {
                    residues.push_back(&res);
                }
            }
        }
    }

    // Test 1: Normal case with water_dist_max = 3.25
    {
        double water_dist_max = 3.25;
        std::vector<const gemmi::Residue*> filtered_residues = 
            coot::filter_residues_by_solvent_contact_gemmi(*res_ref, st, residues, water_dist_max);

        EXPECT_EQ(filtered_residues.size(), 6) 
            << "Expected 6 residues after filter with water_dist_max=3.25";
        
        // Verify all expected residues are present
        std::vector<std::string> expected_residues = {
            "ILE 104", "ILE 99", "LEU 101", "LEU 103", "GLY 98", "PRO 102"
        };
        
        for (const gemmi::Residue* res : filtered_residues) {
            std::string res_name = res->name;
            int seq_num = res->seqid.num.value;
            bool found = false;
            for (const auto& exp : expected_residues) {
                std::string expected_name = exp.substr(0, 3);
                std::string expected_num = exp.substr(4);
                if (res_name == expected_name && std::to_string(seq_num) == expected_num) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << "Residue not found in expected list: " 
                               << res_name << " " << seq_num;
        }
        
        EXPECT_EQ(residues.size() - filtered_residues.size(), 0)
            << "Expected 0 residues to be filtered with water_dist_max=3.25";
    }

    // Test 2: Empty input list
    {
        std::vector<const gemmi::Residue*> empty_residues;
        std::vector<const gemmi::Residue*> result = 
            coot::filter_residues_by_solvent_contact_gemmi(*res_ref, st, empty_residues, 3.25);
        EXPECT_EQ(result.size(), 0) << "Empty input should return empty output";
    }

    // Test 3: Very small water distance (1.5)
    {
        double water_dist_max = 1.5;
        std::vector<const gemmi::Residue*> result = 
            coot::filter_residues_by_solvent_contact_gemmi(*res_ref, st, residues, water_dist_max);
        EXPECT_EQ(result.size(), 6) 
            << "Expected 6 residues after filter with water_dist_max=1.5";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}