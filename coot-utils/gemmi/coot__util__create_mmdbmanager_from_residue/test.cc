#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, CreateMMDBManagerFromResidue) {
    // case: valid residue (ILE 10 in chain A)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        gemmi::CRA found_cra{nullptr, nullptr, nullptr};
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    if (res.name == "ILE" && res.seqid.num.value == 10 && chain.name == "A") {
                        found_cra = {&chain, &res, nullptr};
                        break;
                    }
        
        EXPECT_NE(found_cra.residue, nullptr) << "residue not found";
        EXPECT_NE(found_cra.chain, nullptr) << "chain not found";
        
        if (found_cra.residue && found_cra.chain) {
            gemmi::Structure* mol = coot::util::create_mmdbmanager_from_residue_gemmi(found_cra);
            EXPECT_TRUE(mol != nullptr) << "mol_created should be true";
            
            if (mol) {
                EXPECT_EQ(mol->models.size(), 1);
                
                int nChains = 0;
                int nResidues = 0;
                for (auto& model : mol->models) {
                    nChains += model.chains.size();
                    for (auto& chain : model.chains) {
                        nResidues += chain.residues.size();
                    }
                }
                EXPECT_EQ(nChains, 1);
                EXPECT_EQ(nResidues, 1);
                EXPECT_EQ(found_cra.chain->name, std::string("A"));
                
                delete mol;
            }
        } else {
            EXPECT_TRUE(false) << "residue/chain should not be null";
        }
    }

    // case: null input (no chain)
    {
        gemmi::CRA null_cra{nullptr, nullptr, nullptr};
        gemmi::Structure* mol = coot::util::create_mmdbmanager_from_residue_gemmi(null_cra);
        EXPECT_EQ(mol, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}