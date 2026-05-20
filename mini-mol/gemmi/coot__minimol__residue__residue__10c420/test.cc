#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ResidueConversion) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "PDB loading failed";

    // Case 1: valid residue //A/10
    {
        gemmi::Model& model = st.models[0];
        gemmi::Chain* chain_a = nullptr;
        for (auto& c : model.chains) {
            if (c.name == "A") {
                chain_a = &c;
                break;
            }
        }
        ASSERT_NE(chain_a, nullptr) << "Chain A not found";
        
        gemmi::Residue* res = nullptr;
        for (auto& r : chain_a->residues) {
            if (r.seqid.num.value == 10) {
                res = &r;
                break;
            }
        }
        ASSERT_NE(res, nullptr) << "Residue not found for seqnum 10";
        
        int seqnum = res->seqid.num.value;
        char ins_code = res->seqid.icode;
        std::string resname = res->name;
        
        // Construct minimol residue
        coot::minimol::residue mres(res);
        
        // Check residue name
        EXPECT_EQ(mres.name, "ILE");
        
        // Check sequence number
        EXPECT_EQ(mres.seqnum, 10);
        
        // Check insertion code - gemmi uses ' ' for plain residues
        EXPECT_EQ(mres.ins_code, ' ');
        
        // Check atom count
        EXPECT_EQ(mres.n_atoms(), 8);
        
        // Check first atom
        if (mres.n_atoms() > 0) {
            const coot::minimol::atom &at = mres[0];
            EXPECT_EQ(at.name, " N  ");
            EXPECT_EQ(at.element, " N");
        }
    }

    // Case 2: another residue //A/20
    {
        gemmi::Model& model = st.models[0];
        gemmi::Chain* chain_a = nullptr;
        for (auto& c : model.chains) {
            if (c.name == "A") {
                chain_a = &c;
                break;
            }
        }
        ASSERT_NE(chain_a, nullptr) << "Chain A not found";
        
        gemmi::Residue* res = nullptr;
        for (auto& r : chain_a->residues) {
            if (r.seqid.num.value == 20) {
                res = &r;
                break;
            }
        }
        ASSERT_NE(res, nullptr) << "Residue not found for seqnum 20";
        
        coot::minimol::residue mres(res);
        
        // Check residue name
        EXPECT_EQ(mres.name, "LYS");
        
        // Check atom count
        EXPECT_EQ(mres.n_atoms(), 5);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}