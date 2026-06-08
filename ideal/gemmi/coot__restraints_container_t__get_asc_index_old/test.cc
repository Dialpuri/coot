#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_asc_index_old) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Find chain A, residue 57
    gemmi::Residue* target_res = nullptr;
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 57) {
                    target_res = &res;
                    break;
                }
            }
        }
    }
    ASSERT_NE(target_res, nullptr);

    std::vector<std::pair<bool, gemmi::Residue*>> residues;
    residues.push_back(std::make_pair(false, target_res));

    coot::restraints_container_t_gemmi restraints(residues, &st);

    // Case 1: valid atom - CA in residue 57 chain A
    {
        std::string at_name = "CA";
        int resno = 57;
        std::string chain_id = "A";
        int index = restraints.get_asc_index_old_gemmi(at_name, resno, chain_id);
        EXPECT_EQ(index, 1);
    }

    // Case 2: invalid atom name that doesn't exist in residue 57
    {
        std::string at_name = "ZZZ";
        int resno = 57;
        std::string chain_id = "A";
        int index = restraints.get_asc_index_old_gemmi(at_name, resno, chain_id);
        EXPECT_EQ(index, -1);
    }

    // Case 3: valid atom - CB in residue 57 chain A
    {
        std::string at_name = "CB";
        int resno = 57;
        std::string chain_id = "A";
        int index = restraints.get_asc_index_old_gemmi(at_name, resno, chain_id);
        EXPECT_EQ(index, 4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
