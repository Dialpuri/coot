#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <iostream>
#include "function.hh"

TEST(OracleTest, SpecsToAtomSelection) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // case 1: one residue, atom_mask_mode=0 (all atoms)
    {
        std::vector<coot::residue_spec_t> specs;
        specs.push_back(coot::residue_spec_t("A", 10, ""));

        std::cout << "INPUT residue_specs.size(): " << specs.size() << std::endl;
        std::cout << "INPUT residue_spec[0]: chain=" << specs[0].chain_id
                  << " res_no=" << specs[0].res_no << " ins_code=" << specs[0].ins_code << std::endl;
        std::cout << "INPUT atom_mask_mode: 0" << std::endl;

        int sel_hnd = coot::specs_to_atom_selection_gemmi(specs, st, 0);
        std::cout << "OUTPUT sel_hnd: " << sel_hnd << std::endl;

        EXPECT_GE(sel_hnd, 0);

        if (sel_hnd >= 0) {
            EXPECT_EQ(sel_hnd, 1); // In gemmi we return 1 on success
        }
    }

    // case 2: main-chain mask (atom_mask_mode=1)
    {
        std::vector<coot::residue_spec_t> specs;
        specs.push_back(coot::residue_spec_t("A", 10, ""));

        std::cout << "INPUT residue_specs.size(): " << specs.size() << std::endl;
        std::cout << "INPUT residue_spec[0]: chain=" << specs[0].chain_id
                  << " res_no=" << specs[0].res_no << " ins_code=" << specs[0].ins_code << std::endl;
        std::cout << "INPUT atom_mask_mode: 1 (main-chain)" << std::endl;

        int sel_hnd = coot::specs_to_atom_selection_gemmi(specs, st, 1);
        std::cout << "OUTPUT sel_hnd: " << sel_hnd << std::endl;

        EXPECT_GE(sel_hnd, 0);

        if (sel_hnd >= 0) {
            EXPECT_EQ(sel_hnd, 1);
        }
    }

    // case 3: multiple residues, mask=0
    {
        std::vector<coot::residue_spec_t> specs;
        specs.push_back(coot::residue_spec_t("A", 10, ""));
        specs.push_back(coot::residue_spec_t("A", 11, ""));
        specs.push_back(coot::residue_spec_t("B", 10, ""));

        std::cout << "INPUT residue_specs.size(): " << specs.size() << std::endl;
        for (unsigned int i=0; i<specs.size(); i++)
            std::cout << "INPUT residue_spec[" << i << "]: chain=" << specs[i].chain_id
                      << " res_no=" << specs[i].res_no << " ins_code=" << specs[i].ins_code << std::endl;
        std::cout << "INPUT atom_mask_mode: 0" << std::endl;

        int sel_hnd = coot::specs_to_atom_selection_gemmi(specs, st, 0);
        std::cout << "OUTPUT sel_hnd: " << sel_hnd << std::endl;

        EXPECT_GE(sel_hnd, 0);

        if (sel_hnd >= 0) {
            EXPECT_EQ(sel_hnd, 1);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}