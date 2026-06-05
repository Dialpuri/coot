#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_trp_indole_hydrogens) {
    // case: TRP indole hydrogens on 1f86.pdb
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-dev/coot/reference-structures/1f86.pdb");

        ASSERT_FALSE(st.models.empty());
        gemmi::Model& model = st.models[0];

        ASSERT_TRUE(model.chains.size() > 1);
        gemmi::Chain& chain = model.chains[1];  // GetChain(1) — second chain

        gemmi::Residue* trp_res = nullptr;
        for (gemmi::Residue& res : chain.residues) {
            if (res.name == "TRP") {
                trp_res = &res;
                break;
            }
        }

        ASSERT_NE(trp_res, nullptr);

        std::string cid = "//" + chain.name + "/" + std::to_string(trp_res->seqid.num.value);
        EXPECT_EQ(cid, "//B/41");

        int atoms_before = static_cast<int>(trp_res->atoms.size());
        EXPECT_EQ(atoms_before, 14);

        double bl_HN = 1.01;
        double bl_arom = 1.08;

        EXPECT_NO_THROW(coot::reduce::add_trp_indole_hydrogens_gemmi(*trp_res, bl_HN, bl_arom));

        int atoms_after = static_cast<int>(trp_res->atoms.size());
        EXPECT_EQ(atoms_after, 20);
        EXPECT_EQ(atoms_after - atoms_before, 6);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}