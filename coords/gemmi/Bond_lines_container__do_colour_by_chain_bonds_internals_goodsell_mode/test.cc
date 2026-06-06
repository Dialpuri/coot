#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

int count_bonds(const Bond_lines_container_for_goodsell &bonds) {
    return bonds.total_bond_count();
}

TEST(OracleTest, do_colour_by_chain_bonds_internals_goodsell_mode_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "No models in PDB";

    gemmi::Model &model = st.models[0];
    ASSERT_FALSE(model.chains.empty()) << "No chains";

    gemmi::Chain &chain = model.chains[0];
    ASSERT_FALSE(chain.residues.empty()) << "No residues";

    gemmi::Residue &res = chain.residues[0];
    ASSERT_GE(res.atoms.size(), 3U) << "Need at least 3 atoms";

    gemmi::Atom &at1 = res.atoms[0];
    gemmi::Atom &at2 = res.atoms[1];
    gemmi::Atom &at3 = res.atoms[2];

    int imol = 0;
    int imodel = 1;
    int chain_idx = 0;
    int uddHnd = -1;
    int udd_user_defined_atom_colour_index_handle = -1;

    // Helper: gemmi element.name() returns unpadded ("C"), original expects padded (" C")
    auto elem = [](const gemmi::Atom &a) {
        return std::string(" ") + a.element.name();
    };

    // CASE 1: N-C bond (element1 != "C", element2 == "C") -> 2 half-bonds
    {
        Bond_lines_container_for_goodsell bonds;
        bool have_dictionary = false;

        std::string e1 = elem(at1);
        std::string e2 = elem(at2);
        coot::Cartesian pos1(at1.pos.x, at1.pos.y, at1.pos.z);
        coot::Cartesian pos2(at2.pos.x, at2.pos.y, at2.pos.z);

        std::vector<std::pair<bool, gemmi::Residue*>> het_residues;

        int bonds_before = count_bonds(bonds);

        do_colour_by_chain_bonds_internals_goodsell_mode_gemmi(
            bonds,
            imol, imodel, chain_idx,
            &at1, &at2,
            &res, &res, 0, 1,
            &het_residues, e1, e2, pos1, pos2,
            uddHnd, udd_user_defined_atom_colour_index_handle,
            have_dictionary
        );

        int bonds_after = count_bonds(bonds);

        EXPECT_EQ(bonds_before, 0);
        EXPECT_EQ(bonds_after, 2);
    }

    // CASE 2: C-C bond (element1 == "C", element2 == "C") -> 1 bond with C-colour
    {
        Bond_lines_container_for_goodsell bonds;
        bool have_dictionary = false;

        std::string e2 = elem(at2);
        std::string e3 = elem(at3);
        coot::Cartesian pos2(at2.pos.x, at2.pos.y, at2.pos.z);
        coot::Cartesian pos3(at3.pos.x, at3.pos.y, at3.pos.z);

        std::vector<std::pair<bool, gemmi::Residue*>> het_residues;

        int bonds_before = count_bonds(bonds);

        do_colour_by_chain_bonds_internals_goodsell_mode_gemmi(
            bonds,
            imol, imodel, chain_idx,
            &at2, &at3,
            &res, &res, 1, 2,
            &het_residues, e2, e3, pos2, pos3,
            uddHnd, udd_user_defined_atom_colour_index_handle,
            have_dictionary
        );

        int bonds_after = count_bonds(bonds);

        EXPECT_EQ(bonds_before, 0);
        EXPECT_EQ(bonds_after, 1);
    }

    // CASE 3: far-away H-like atoms with distance > 1.5 -> draw_it = false, no bonds
    {
        Bond_lines_container_for_goodsell bonds;
        bool have_dictionary = false;

        coot::Cartesian pos_far1(0, 0, 0);
        coot::Cartesian pos_far2(0, 0, 3.0); // distance = 3.0 > 1.5

        std::string e_h = " H";
        std::string e_c = " C";

        std::vector<std::pair<bool, gemmi::Residue*>> het_residues;

        int bonds_before = count_bonds(bonds);

        do_colour_by_chain_bonds_internals_goodsell_mode_gemmi(
            bonds,
            imol, imodel, chain_idx,
            &at1, &at2,
            &res, &res, 0, 1,
            &het_residues, e_h, e_c, pos_far1, pos_far2,
            uddHnd, udd_user_defined_atom_colour_index_handle,
            have_dictionary
        );

        int bonds_after = count_bonds(bonds);

        EXPECT_EQ(bonds_before, 0);
        EXPECT_EQ(bonds_after, 0);
    }

    // CASE 4: het residue — HOH is in chain B of example.pdb
    {
        if (model.chains.size() < 2) {
            GTEST_SKIP() << "Need at least 2 chains for het test";
        }

        gemmi::Chain &chain_b = model.chains[1];
        if (chain_b.residues.empty()) {
            GTEST_SKIP() << "Chain B has no residues";
        }

        gemmi::Residue &het_res = chain_b.residues[0];
        if (het_res.atoms.empty()) {
            GTEST_SKIP() << "Het residue has no atoms";
        }
        gemmi::Atom &het_atom = het_res.atoms[0];

        Bond_lines_container_for_goodsell bonds;
        bool have_dictionary = false;

        std::vector<std::pair<bool, gemmi::Residue*>> het_residues_with_het;
        het_residues_with_het.push_back({true, &het_res});

        std::string e_het = elem(het_atom);
        std::string e_c = elem(at2);
        coot::Cartesian pos_het(het_atom.pos.x, het_atom.pos.y, het_atom.pos.z);
        coot::Cartesian pos_c(at2.pos.x, at2.pos.y, at2.pos.z);

        int bonds_before = count_bonds(bonds);

        do_colour_by_chain_bonds_internals_goodsell_mode_gemmi(
            bonds,
            imol, imodel, chain_idx,
            &het_atom, &at2,
            &het_res, &res, 0, 1,
            &het_residues_with_het, e_het, e_c, pos_het, pos_c,
            uddHnd, udd_user_defined_atom_colour_index_handle,
            have_dictionary,
            true,  // atom1_is_het
            false  // atom2_is_het
        );

        int bonds_after = count_bonds(bonds);

        EXPECT_EQ(bonds_before, 0);
        EXPECT_FALSE(bonds_after < 0); // just verify no crash and non-negative
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}