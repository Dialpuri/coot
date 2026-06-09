#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_his_ring_C_Hs) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Case 1: valid HIS residue (//A/60 is HIS)
    {
        gemmi::Residue* his_res = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 60) {
                    his_res = &res;
                    break;
                }
            }
        }
        ASSERT_NE(his_res, nullptr);

        double bl_arom = 1.08;

        int atoms_before = static_cast<int>(his_res->atoms.size());
        std::string chain_id = "A";
        int seqnum = his_res->seqid.num.value;
        std::string resname = his_res->name;

        EXPECT_EQ(chain_id, "A");
        EXPECT_EQ(seqnum, 60);
        EXPECT_EQ(resname, "HIS");
        EXPECT_NEAR(bl_arom, 1.08, 1e-4);
        EXPECT_EQ(atoms_before, 10);

        coot::reduce::add_his_ring_C_Hs_gemmi(*his_res, bl_arom);

        int atoms_after = static_cast<int>(his_res->atoms.size());
        EXPECT_EQ(atoms_after, 12);
        EXPECT_EQ((atoms_after - atoms_before), 2);
    }

    // Case 2: non-HIS residue (MET at //A/1)
    {
        gemmi::Residue* non_his = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 1) {
                    non_his = &res;
                    break;
                }
            }
        }
        ASSERT_NE(non_his, nullptr);

        double bl_arom = 1.08;

        int atoms_before2 = static_cast<int>(non_his->atoms.size());
        std::string chain_id2 = "A";
        int seqnum2 = non_his->seqid.num.value;
        std::string resname2 = non_his->name;

        EXPECT_EQ(chain_id2, "A");
        EXPECT_EQ(seqnum2, 1);
        EXPECT_EQ(resname2, "MET");
        EXPECT_NEAR(bl_arom, 1.08, 1e-4);
        EXPECT_EQ(atoms_before2, 8);

        coot::reduce::add_his_ring_C_Hs_gemmi(*non_his, bl_arom);

        int atoms_after2 = static_cast<int>(non_his->atoms.size());
        EXPECT_EQ(atoms_after2, 8);
        EXPECT_EQ((atoms_after2 - atoms_before2), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
