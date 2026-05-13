#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, make_flat_ligand_name_map) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];
    gemmi::Residue* res = nullptr;
    
    for (auto& r : chain.residues) {
        if (r.seqid.num.value == 10) {
            res = &r;
            break;
        }
    }
    
    ASSERT_NE(res, nullptr) << "Residue 10 not found";
    
    std::map<std::string, std::string> name_map = pli::make_flat_ligand_name_map_gemmi(*res);
    
    EXPECT_EQ(name_map.size(), 0u);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}