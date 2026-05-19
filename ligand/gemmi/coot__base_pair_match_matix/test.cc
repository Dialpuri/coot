#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, BasePairMatchMatrix) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case: valid guanine residues with identical atoms
    {
        gemmi::Model& model = st.models[0];
        
        // Find residues A/10 and A/20
        gemmi::Residue *res1 = nullptr;
        gemmi::Residue *res2 = nullptr;
        
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 10) res1 = &res;
                    if (res.seqid.num.value == 20) res2 = &res;
                }
            }
        }
        
        ASSERT_NE(res1, nullptr);
        ASSERT_NE(res2, nullptr);
        
        res1->name = "G";
        res2->name = "G";
        
        // Add base atoms to res1
        gemmi::Atom at1;
        at1.name = " N9 ";
        at1.pos = gemmi::Position(10.0, 10.0, 10.0);
        at1.occ = 1.0;
        at1.b_iso = 20.0;
        res1->atoms.push_back(at1);
        
        gemmi::Atom at2;
        at2.name = " C8 ";
        at2.pos = gemmi::Position(11.0, 10.0, 10.0);
        at2.occ = 1.0;
        at2.b_iso = 20.0;
        res1->atoms.push_back(at2);
        
        gemmi::Atom at3;
        at3.name = " N7 ";
        at3.pos = gemmi::Position(11.0, 11.0, 10.0);
        at3.occ = 1.0;
        at3.b_iso = 20.0;
        res1->atoms.push_back(at3);
        
        gemmi::Atom at4;
        at4.name = " C5 ";
        at4.pos = gemmi::Position(10.0, 11.0, 10.0);
        at4.occ = 1.0;
        at4.b_iso = 20.0;
        res1->atoms.push_back(at4);
        
        gemmi::Atom at5;
        at5.name = " C4 ";
        at5.pos = gemmi::Position(10.0, 10.5, 11.0);
        at5.occ = 1.0;
        at5.b_iso = 20.0;
        res1->atoms.push_back(at5);
        
        gemmi::Atom at6;
        at6.name = " N1 ";
        at6.pos = gemmi::Position(11.0, 10.5, 11.0);
        at6.occ = 1.0;
        at6.b_iso = 20.0;
        res1->atoms.push_back(at6);
        
        gemmi::Atom at7;
        at7.name = " C2 ";
        at7.pos = gemmi::Position(11.0, 11.5, 11.0);
        at7.occ = 1.0;
        at7.b_iso = 20.0;
        res1->atoms.push_back(at7);
        
        gemmi::Atom at8;
        at8.name = " N3 ";
        at8.pos = gemmi::Position(10.0, 11.5, 11.0);
        at8.occ = 1.0;
        at8.b_iso = 20.0;
        res1->atoms.push_back(at8);
        
        gemmi::Atom at9;
        at9.name = " C6 ";
        at9.pos = gemmi::Position(10.5, 10.0, 10.0);
        at9.occ = 1.0;
        at9.b_iso = 20.0;
        res1->atoms.push_back(at9);
        
        gemmi::Atom at10;
        at10.name = " N6 ";
        at10.pos = gemmi::Position(10.5, 10.0, 9.0);
        at10.occ = 1.0;
        at10.b_iso = 20.0;
        res1->atoms.push_back(at10);
        
        // Add the same atoms to res2 (identical position - should give identity transform)
        gemmi::Atom at1b;
        at1b.name = " N9 ";
        at1b.pos = gemmi::Position(10.0, 10.0, 10.0);
        at1b.occ = 1.0;
        at1b.b_iso = 20.0;
        res2->atoms.push_back(at1b);
        
        gemmi::Atom at2b;
        at2b.name = " C8 ";
        at2b.pos = gemmi::Position(11.0, 10.0, 10.0);
        at2b.occ = 1.0;
        at2b.b_iso = 20.0;
        res2->atoms.push_back(at2b);
        
        gemmi::Atom at3b;
        at3b.name = " N7 ";
        at3b.pos = gemmi::Position(11.0, 11.0, 10.0);
        at3b.occ = 1.0;
        at3b.b_iso = 20.0;
        res2->atoms.push_back(at3b);
        
        gemmi::Atom at4b;
        at4b.name = " C5 ";
        at4b.pos = gemmi::Position(10.0, 11.0, 10.0);
        at4b.occ = 1.0;
        at4b.b_iso = 20.0;
        res2->atoms.push_back(at4b);
        
        gemmi::Atom at5b;
        at5b.name = " C4 ";
        at5b.pos = gemmi::Position(10.0, 10.5, 11.0);
        at5b.occ = 1.0;
        at5b.b_iso = 20.0;
        res2->atoms.push_back(at5b);
        
        gemmi::Atom at6b;
        at6b.name = " N1 ";
        at6b.pos = gemmi::Position(11.0, 10.5, 11.0);
        at6b.occ = 1.0;
        at6b.b_iso = 20.0;
        res2->atoms.push_back(at6b);
        
        gemmi::Atom at7b;
        at7b.name = " C2 ";
        at7b.pos = gemmi::Position(11.0, 11.5, 11.0);
        at7b.occ = 1.0;
        at7b.b_iso = 20.0;
        res2->atoms.push_back(at7b);
        
        gemmi::Atom at8b;
        at8b.name = " N3 ";
        at8b.pos = gemmi::Position(10.0, 11.5, 11.0);
        at8b.occ = 1.0;
        at8b.b_iso = 20.0;
        res2->atoms.push_back(at8b);
        
        gemmi::Atom at9b;
        at9b.name = " C6 ";
        at9b.pos = gemmi::Position(10.5, 10.0, 10.0);
        at9b.occ = 1.0;
        at9b.b_iso = 20.0;
        res2->atoms.push_back(at9b);
        
        gemmi::Atom at10b;
        at10b.name = " N6 ";
        at10b.pos = gemmi::Position(10.5, 10.0, 9.0);
        at10b.occ = 1.0;
        at10b.b_iso = 20.0;
        res2->atoms.push_back(at10b);
        
        std::pair<bool, clipper::RTop_orth> result = coot::base_pair_match_matix_gemmi(res1, res2);
        
        EXPECT_TRUE(result.first);
        
        if (result.first) {
            clipper::Mat33<double> rot = result.second.rot();
            clipper::Vec3<double> trans = result.second.trn();
            EXPECT_NEAR(trans[0], 0.0, 1e-4);
            EXPECT_NEAR(trans[1], 0.0, 1e-4);
            EXPECT_NEAR(trans[2], 0.0, 1e-4);
        }
    }
    
    // case: invalid residue names
    {
        gemmi::Model& model = st.models[0];
        
        gemmi::Residue *res1 = nullptr;
        gemmi::Residue *res2 = nullptr;
        
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 10) res1 = &res;
                    if (res.seqid.num.value == 20) res2 = &res;
                }
            }
        }
        
        ASSERT_NE(res1, nullptr);
        ASSERT_NE(res2, nullptr);
        
        res1->name = "XXX";
        res2->name = "YYY";
        
        std::pair<bool, clipper::RTop_orth> result = coot::base_pair_match_matix_gemmi(res1, res2);
        
        EXPECT_FALSE(result.first);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}