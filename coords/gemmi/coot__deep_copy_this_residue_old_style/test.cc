#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, DeepCopyResidue) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model &model = st.models[0];
    gemmi::Chain &chain = model.chains[0];  // chain A
    gemmi::Residue &res = chain.residues[9];  // residue 10 (0-indexed)

    // Build CRA for the residue we want to copy
    gemmi::CRA cra{&chain, &res, nullptr};

    // Case 1: whole_residue_flag=1 (copy all atoms regardless of altLoc), embed_in_chain_flag=true
    {
        std::cout << "INPUT residue: " << cra.chain->name << " " << cra.residue->seqid.num.value << " " << cra.residue->seqid.icode << std::endl;
        std::cout << "INPUT altconf: \"\"" << std::endl;
        std::cout << "INPUT whole_residue_flag: 1" << std::endl;
        std::cout << "INPUT embed_in_chain_flag: 1" << std::endl;

        int n_atoms_before = (int)cra.residue->atoms.size();
        std::cout << "INPUT n_atoms_before: " << n_atoms_before << std::endl;

        gemmi::Residue *rres = coot::deep_copy_this_residue_old_style_gemmi(cra, "", 1, 1);

        std::cout << "OUTPUT rres: " << (rres ? "non-null" : "null") << std::endl;
        EXPECT_NE(rres, nullptr);
        if (rres) {
            std::cout << "OUTPUT copied_n_atoms: " << (int)rres->atoms.size() << std::endl;
            EXPECT_EQ((int)rres->atoms.size(), 8);
            delete rres;
        }
    }

    // Case 2: whole_residue_flag=0, altconf="" (copy atoms with altLoc ' ' or matching altconf)
    {
        std::cout << "INPUT residue (case2): " << cra.chain->name << " " << cra.residue->seqid.num.value << std::endl;
        std::cout << "INPUT altconf (case2): \"\"" << std::endl;
        std::cout << "INPUT whole_residue_flag (case2): 0" << std::endl;
        std::cout << "INPUT embed_in_chain_flag (case2): 0" << std::endl;

        gemmi::Residue *rres2 = coot::deep_copy_this_residue_old_style_gemmi(cra, "", 0, 0);
        std::cout << "OUTPUT rres2: " << (rres2 ? "non-null" : "null") << std::endl;
        EXPECT_NE(rres2, nullptr);
        if (rres2) {
            std::cout << "OUTPUT copied_n_atoms2: " << (int)rres2->atoms.size() << std::endl;
            EXPECT_EQ((int)rres2->atoms.size(), 8);
            delete rres2;
        }
    }

    // Case 3: residue with no atoms (edge case — should return NULL)
    {
        std::cout << "INPUT residue (empty): A 999 \"\" " << std::endl;
        gemmi::Residue empty_res;
        empty_res.name = "ALA";
        empty_res.seqid.num.value = 999;
        empty_res.seqid.icode = ' ';
        gemmi::CRA empty_cra{nullptr, &empty_res, nullptr};
        gemmi::Residue *rres3 = coot::deep_copy_this_residue_old_style_gemmi(empty_cra, "", 1, 0);
        std::cout << "OUTPUT rres_empty: " << (rres3 ? "non-null" : "null") << std::endl;
        EXPECT_EQ(rres3, nullptr);
        if (rres3) delete rres3;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}