#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cifdoc.hpp>
#include <gemmi/read_cif.hpp>
#include <gemmi/unitcell.hpp>
#include <gemmi/math.hpp>
#include <clipper/core/clipper_types.h>
#include <clipper/core/spacegroup.h>
#include "function.hh"
#include <fstream>
#include <string>
#include <vector>

TEST(OracleTest, read_coordinates) {
    // ========================================
    // Case 1: CIF with atoms
    // ========================================
    {
        std::string cif_content =
            "data_test\n"
            "_cell_length_a   10.0\n"
            "_cell_length_b   10.0\n"
            "_cell_length_c   10.0\n"
            "_cell_angle_alpha 90.0\n"
            "_cell_angle_beta  90.0\n"
            "_cell_angle_gamma 90.0\n"
            "_symmetry_space_group_name_H-M 'P 1'\n"
            "loop_\n"
            "_symmetry_equiv_pos_as_xyz\n"
            "'x, y, z'\n"
            "loop_\n"
            "_atom_site_label\n"
            "_atom_site_fract_x\n"
            "_atom_site_fract_y\n"
            "_atom_site_fract_z\n"
            "_atom_site_type_symbol\n"
            "_atom_site_occupancy\n"
            "_atom_site_U_iso_or_equiv\n"
            "C1  0.1 0.2 0.3 C 1.0 0.05\n"
            "C2  0.2 0.3 0.4 C 1.0 0.06\n"
            "O1  0.3 0.4 0.5 O 1.0 0.07\n"
            ;

        std::string cif_path = "/tmp/test_sm_cif.cif";
        {
            std::ofstream ofs(cif_path);
            ofs << cif_content;
        }

        // Read CIF with gemmi
        gemmi::cif::Document doc = gemmi::read_cif_gz(cif_path);
        EXPECT_FALSE(doc.blocks.empty());

        gemmi::cif::Block& bloc = doc.blocks[0];

        // Get cell from CIF
        const std::string* a_ptr = bloc.find_value("_cell_length_a");
        const std::string* b_ptr = bloc.find_value("_cell_length_b");
        const std::string* c_ptr = bloc.find_value("_cell_length_c");
        const std::string* alpha_ptr = bloc.find_value("_cell_angle_alpha");
        const std::string* beta_ptr = bloc.find_value("_cell_angle_beta");
        const std::string* gamma_ptr = bloc.find_value("_cell_angle_gamma");

        double a = a_ptr ? std::stod(*a_ptr) : 10.0;
        double b = b_ptr ? std::stod(*b_ptr) : 10.0;
        double c = c_ptr ? std::stod(*c_ptr) : 10.0;
        double al = alpha_ptr ? std::stod(*alpha_ptr) : 90.0;
        double be = beta_ptr ? std::stod(*beta_ptr) : 90.0;
        double ga = gamma_ptr ? std::stod(*gamma_ptr) : 90.0;

        clipper::Cell_descr cell_descr(a, b, c,
                                       clipper::Util::d2rad(al),
                                       clipper::Util::d2rad(be),
                                       clipper::Util::d2rad(ga));
        clipper::Cell cell(cell_descr);

        // Verify cell format matches expected
        EXPECT_EQ(std::string(cell.format()), " Cell (    10,    10,    10,    90,    90,    90)");

        // Get space group from symmetry operations
        gemmi::cif::Column sym_col = bloc.find_loop("_symmetry_equiv_pos_as_xyz");
        std::vector<std::string> symm_strings;
        if (sym_col.length() > 0) {
            for (int i = 0; i < sym_col.length(); ++i) {
                symm_strings.push_back(sym_col.at(i));
            }
        }
        EXPECT_EQ(symm_strings.size(), 1u);

        // Build space group from symmetry operations
        std::string symmetry_ops;
        for (const auto& s : symm_strings) {
            symmetry_ops += s;
            symmetry_ops += " ; ";
        }
        clipper::Spgr_descr spg_descr(symmetry_ops, clipper::Spgr_descr::Symops);
        clipper::Spacegroup spg;
        spg.init(spg_descr);
        EXPECT_EQ(std::string(spg.symbol_xhm()), "P 1");

        // Call read_coordinates
        std::vector<gemmi::Atom> atoms = coot::smcif::read_coordinates_gemmi(cif_path, cell, spg);
        EXPECT_EQ(atoms.size(), 3u);

        // Check atom 0
        EXPECT_NE(&atoms[0], nullptr);
        EXPECT_EQ(atoms[0].name, "C1");
        EXPECT_EQ(std::string(atoms[0].element.name()), "C");

        // Check atom 1
        EXPECT_NE(&atoms[1], nullptr);
        EXPECT_EQ(atoms[1].name, "C2");
        EXPECT_EQ(std::string(atoms[1].element.name()), "C");

        // Check atom 2
        EXPECT_NE(&atoms[2], nullptr);
        EXPECT_EQ(atoms[2].name, "O1");
        EXPECT_EQ(std::string(atoms[2].element.name()), "O");
    }

    // ========================================
    // Case 2: CIF with no atom loop (edge case)
    // ========================================
    {
        std::string cif_empty =
            "data_test\n"
            "_cell_length_a   10.0\n"
            "_cell_length_b   10.0\n"
            "_cell_length_c   10.0\n"
            "_cell_angle_alpha 90.0\n"
            "_cell_angle_beta  90.0\n"
            "_cell_angle_gamma 90.0\n"
            "_symmetry_space_group_name_H-M 'P 1'\n"
            "loop_\n"
            "_symmetry_equiv_pos_as_xyz\n"
            "'x, y, z'\n"
            ;

        std::string cif_empty_path = "/tmp/test_sm_cif_empty.cif";
        {
            std::ofstream ofs(cif_empty_path);
            ofs << cif_empty;
        }

        // Read CIF with gemmi
        gemmi::cif::Document doc = gemmi::read_cif_gz(cif_empty_path);
        EXPECT_FALSE(doc.blocks.empty());

        gemmi::cif::Block& bloc = doc.blocks[0];

        // Get cell from CIF
        const std::string* a_ptr = bloc.find_value("_cell_length_a");
        const std::string* b_ptr = bloc.find_value("_cell_length_b");
        const std::string* c_ptr = bloc.find_value("_cell_length_c");
        const std::string* alpha_ptr = bloc.find_value("_cell_angle_alpha");
        const std::string* beta_ptr = bloc.find_value("_cell_angle_beta");
        const std::string* gamma_ptr = bloc.find_value("_cell_angle_gamma");

        double a = a_ptr ? std::stod(*a_ptr) : 10.0;
        double b = b_ptr ? std::stod(*b_ptr) : 10.0;
        double c = c_ptr ? std::stod(*c_ptr) : 10.0;
        double al = alpha_ptr ? std::stod(*alpha_ptr) : 90.0;
        double be = beta_ptr ? std::stod(*beta_ptr) : 90.0;
        double ga = gamma_ptr ? std::stod(*gamma_ptr) : 90.0;

        clipper::Cell_descr cell_descr(a, b, c,
                                       clipper::Util::d2rad(al),
                                       clipper::Util::d2rad(be),
                                       clipper::Util::d2rad(ga));
        clipper::Cell cell(cell_descr);

        EXPECT_EQ(std::string(cell.format()), " Cell (    10,    10,    10,    90,    90,    90)");

        // Get space group
        gemmi::cif::Column sym_col = bloc.find_loop("_symmetry_equiv_pos_as_xyz");
        std::vector<std::string> symm_strings;
        if (sym_col.length() > 0) {
            for (int i = 0; i < sym_col.length(); ++i) {
                symm_strings.push_back(sym_col.at(i));
            }
        }
        EXPECT_EQ(symm_strings.size(), 1u);

        std::string symmetry_ops;
        for (const auto& s : symm_strings) {
            symmetry_ops += s;
            symmetry_ops += " ; ";
        }
        clipper::Spgr_descr spg_descr(symmetry_ops, clipper::Spgr_descr::Symops);
        clipper::Spacegroup spg;
        spg.init(spg_descr);
        EXPECT_EQ(std::string(spg.symbol_xhm()), "P 1");

        // Call read_coordinates
        std::vector<gemmi::Atom> atoms = coot::smcif::read_coordinates_gemmi(cif_empty_path, cell, spg);
        EXPECT_EQ(atoms.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}