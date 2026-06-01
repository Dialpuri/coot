#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static int count_atoms_manual(const gemmi::Structure& st) {
    int total = 0;
    for (const gemmi::Model& model : st.models)
        for (const gemmi::Chain& chain : model.chains)
            for (const gemmi::Residue& res : chain.residues)
                total += static_cast<int>(res.atoms.size());
    return total;
}

TEST(OracleTest, sort_residues) {
    // Case 1: sort_residues on a valid molecule
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        EXPECT_EQ(st.models.size(), 1);

        gemmi::Model &model = st.models[0];
        int nchains = model.chains.size();

        EXPECT_EQ(st.models.size(), 1);
        EXPECT_EQ(nchains, 2);

        // Capture first 5 residue seqnums from chain 0 BEFORE
        gemmi::Chain &chain = model.chains[0];
        std::string chain_id = chain.name;
        int nres = chain.residues.size();

        EXPECT_EQ(chain_id, "A");
        EXPECT_EQ(nres, 267);

        int nres_cap = std::min(nres, 5);
        std::vector<int> before_seqnums;
        for (int i = 0; i < nres_cap; i++) {
            gemmi::Residue &res = chain.residues[i];
            before_seqnums.push_back(res.seqid.num.value);
        }
        // BEFORE first_5_residue_seqnums: 1,2,3,4,5
        EXPECT_EQ(before_seqnums.size(), 5);
        EXPECT_EQ(before_seqnums[0], 1);
        EXPECT_EQ(before_seqnums[1], 2);
        EXPECT_EQ(before_seqnums[2], 3);
        EXPECT_EQ(before_seqnums[3], 4);
        EXPECT_EQ(before_seqnums[4], 5);

        int total_atoms_before = count_atoms_manual(st);
        EXPECT_EQ(total_atoms_before, 2107);

        // Call sort_residues
        EXPECT_NO_THROW(coot::sort_residues_gemmi(&st));

        // Capture first 5 residue seqnums from chain 0 AFTER
        chain = model.chains[0];
        nres = chain.residues.size();
        nres_cap = std::min(nres, 5);
        std::vector<int> after_seqnums;
        for (int i = 0; i < nres_cap; i++) {
            gemmi::Residue &res = chain.residues[i];
            after_seqnums.push_back(res.seqid.num.value);
        }
        // AFTER first_5_residue_seqnums: 1,2,3,4,5
        EXPECT_EQ(after_seqnums.size(), 5);
        EXPECT_EQ(after_seqnums[0], 1);
        EXPECT_EQ(after_seqnums[1], 2);
        EXPECT_EQ(after_seqnums[2], 3);
        EXPECT_EQ(after_seqnums[3], 4);
        EXPECT_EQ(after_seqnums[4], 5);

        int total_atoms_after = count_atoms_manual(st);
        EXPECT_EQ(total_atoms_after, 2107);
    }

    // Case 2: null pointer — triggers the null-molecule guard
    {
        EXPECT_NO_THROW(coot::sort_residues_gemmi(nullptr));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}