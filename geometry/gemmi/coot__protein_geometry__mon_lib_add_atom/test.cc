#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, mon_lib_add_atom) {
    int imol_enc = 0;
    int read_number = 0;
    std::vector<std::pair<int, coot::dictionary_residue_restraints_t>> dict_res_restraints;

    // Record initial state
    size_t initial_size = dict_res_restraints.size();
    EXPECT_EQ(initial_size, size_t(0));

    // --- Case 1: Add first atom for a new residue (ALA, atom N) ---
    {
        std::string comp_id = "ALA";
        std::string atom_id = "N";
        std::string atom_id_4c = "N   ";
        std::string type_symbol = "N";
        std::string type_energy = "N";
        std::pair<bool, double> partial_charge(false, 0.0);
        std::pair<bool, int> formal_charge(false, 0);
        coot::dict_atom::aromaticity_t arom = coot::dict_atom::NON_AROMATIC;
        std::pair<bool, clipper::Coord_orth> model_pos(false, clipper::Coord_orth(0,0,0));
        std::pair<bool, clipper::Coord_orth> model_pos_ideal(true, clipper::Coord_orth(0.0, 0.0, 1.460));

        size_t before = dict_res_restraints.size();
        EXPECT_EQ(before, size_t(0));

        coot::mon_lib_add_atom_gemmi(dict_res_restraints, read_number, comp_id, imol_enc, atom_id, atom_id_4c, type_symbol,
                              type_energy, partial_charge, formal_charge, arom,
                              model_pos, model_pos_ideal);

        size_t after = dict_res_restraints.size();
        EXPECT_EQ(after, size_t(1));
        EXPECT_EQ(dict_res_restraints[after - 1].second.atom_info.size(), size_t(1));
    }

    // --- Case 2: Add second atom to same residue (ALA, atom CA) ---
    {
        std::string comp_id = "ALA";
        std::string atom_id = "CA";
        std::string atom_id_4c = "CA  ";
        std::string type_symbol = "C";
        std::string type_energy = "C";
        std::pair<bool, double> partial_charge(false, 0.0);
        std::pair<bool, int> formal_charge(false, 0);
        coot::dict_atom::aromaticity_t arom = coot::dict_atom::NON_AROMATIC;
        std::pair<bool, clipper::Coord_orth> model_pos(true, clipper::Coord_orth(1.458, -0.0, 0.0));
        std::pair<bool, clipper::Coord_orth> model_pos_ideal(true, clipper::Coord_orth(2.265, -0.0, 0.0));

        size_t before = dict_res_restraints.size();
        EXPECT_EQ(before, size_t(1));

        coot::mon_lib_add_atom_gemmi(dict_res_restraints, read_number, comp_id, imol_enc, atom_id, atom_id_4c, type_symbol,
                              type_energy, partial_charge, formal_charge, arom,
                              model_pos, model_pos_ideal);

        size_t after = dict_res_restraints.size();
        EXPECT_EQ(after, size_t(1));
        EXPECT_EQ(dict_res_restraints[after - 1].second.atom_info.size(), size_t(2));
    }

    // --- Case 3: Add atom for a different residue (GLY, atom N) ---
    {
        std::string comp_id = "GLY";
        std::string atom_id = "N";
        std::string atom_id_4c = "N   ";
        std::string type_symbol = "N";
        std::string type_energy = "N";
        std::pair<bool, double> partial_charge(true, -0.35);
        std::pair<bool, int> formal_charge(true, 0);
        coot::dict_atom::aromaticity_t arom = coot::dict_atom::UNASSIGNED;
        std::pair<bool, clipper::Coord_orth> model_pos(false, clipper::Coord_orth(0,0,0));
        std::pair<bool, clipper::Coord_orth> model_pos_ideal(false, clipper::Coord_orth(0,0,0));

        size_t before = dict_res_restraints.size();
        EXPECT_EQ(before, size_t(1));

        coot::mon_lib_add_atom_gemmi(dict_res_restraints, read_number, comp_id, imol_enc, atom_id, atom_id_4c, type_symbol,
                              type_energy, partial_charge, formal_charge, arom,
                              model_pos, model_pos_ideal);

        size_t after = dict_res_restraints.size();
        EXPECT_EQ(after, size_t(2));
        EXPECT_EQ(dict_res_restraints[after - 1].second.atom_info.size(), size_t(1));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
