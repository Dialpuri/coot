#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, find_glycosidic_linkage_type_with_order_switch) {
    // Load the protein-ligand complex
    std::string cif_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif";
    gemmi::Structure st = gemmi::read_structure_file(cif_path);

    // Collect ASN residues with their parent chain
    std::vector<std::pair<gemmi::Chain*, gemmi::Residue*>> asn_pairs;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& residue : chain.residues) {
                if (residue.name == "ASN") {
                    asn_pairs.push_back({&chain, &residue});
                }
            }
        }
    }

    // Case 3: Non-glycosylated pair (ASN-ASN)
    // Oracle output: asn1: A/3/ASN, asn2: A/23/ASN, linkage_type: "", order_switch: false
    if (asn_pairs.size() >= 2) {
        // Create CRA objects for the two ASN residues
        gemmi::Chain* chain1 = asn_pairs[0].first;
        gemmi::Residue* res1 = asn_pairs[0].second;
        gemmi::Chain* chain2 = asn_pairs[1].first;
        gemmi::Residue* res2 = asn_pairs[1].second;
        
        gemmi::CRA cra1{chain1, res1, nullptr};
        gemmi::CRA cra2{chain2, res2, nullptr};
        
        std::pair<std::string, bool> result = 
            coot::protein_geometry::find_glycosidic_linkage_type_with_order_switch_gemmi(cra1, cra2);
        EXPECT_EQ(result.first, "") << "linkage_type should be empty for ASN-ASN pair";
        EXPECT_EQ(result.second, false) << "order_switch should be false for ASN-ASN pair";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}