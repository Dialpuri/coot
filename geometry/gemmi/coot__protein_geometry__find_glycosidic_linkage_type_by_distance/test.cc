#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, FindGlycosidicLinkageTypeByDistance) {
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];
        
        gemmi::CRA cra1, cra2;
        for (auto& chain : model.chains)
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 10 && res.name == "ILE") cra1 = {&chain, &res, nullptr};
                if (res.seqid.num.value == 20 && res.name == "LYS") cra2 = {&chain, &res, nullptr};
            }
        
        std::string result = coot::protein_geometry::find_glycosidic_linkage_type_by_distance_gemmi(cra1, cra2);
        EXPECT_EQ(result, "");
    }
    
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];
        
        gemmi::CRA cra1, cra2;
        for (auto& chain : model.chains)
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 10 && res.name == "ILE") cra1 = {&chain, &res, nullptr};
                if (res.seqid.num.value == 11 && res.name == "GLY") cra2 = {&chain, &res, nullptr};
            }
        
        std::string result = coot::protein_geometry::find_glycosidic_linkage_type_by_distance_gemmi(cra1, cra2);
        EXPECT_EQ(result, "");
    }
    
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];
        
        gemmi::CRA cra1, cra2;
        for (auto& chain : model.chains)
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 11 && res.name == "GLY") cra1 = {&chain, &res, nullptr};
                if (res.seqid.num.value == 10 && res.name == "ILE") cra2 = {&chain, &res, nullptr};
            }
        
        std::string result = coot::protein_geometry::find_glycosidic_linkage_type_by_distance_gemmi(cra1, cra2);
        EXPECT_EQ(result, "");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}