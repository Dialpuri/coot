#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <gemmi/cifdoc.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__protein_geometry__link_add_angle/gemmi/function.hh"

namespace coot {
namespace protein_geometry {

void link_angle_gemmi(gemmi::cif::Loop& loop,
                      std::vector<dictionary_residue_link_restraints_t>& dict_link_res_restraints) {
    for (size_t j = 0; j < loop.length(); j++) {
        std::string link_id;
        std::string atom_id_1, atom_id_2, atom_id_3;
        double value_angle, value_angle_esd;
        int atom_1_comp_id, atom_2_comp_id, atom_3_comp_id;

        int ierr_tot = 0;

        int link_id_col = loop.find_tag("link_id");
        if (link_id_col >= 0) {
            link_id = gemmi::cif::as_string(loop.val(j, link_id_col));
        } else {
            ierr_tot++;
        }

        int atom_id_1_col = loop.find_tag("atom_id_1");
        if (atom_id_1_col >= 0) {
            atom_id_1 = gemmi::cif::as_string(loop.val(j, atom_id_1_col));
        } else {
            ierr_tot++;
        }

        int atom_id_2_col = loop.find_tag("atom_id_2");
        if (atom_id_2_col >= 0) {
            atom_id_2 = gemmi::cif::as_string(loop.val(j, atom_id_2_col));
        } else {
            ierr_tot++;
        }

        int atom_id_3_col = loop.find_tag("atom_id_3");
        if (atom_id_3_col >= 0) {
            atom_id_3 = gemmi::cif::as_string(loop.val(j, atom_id_3_col));
        } else {
            ierr_tot++;
        }

        int col1 = loop.find_tag("atom_1_comp_id");
        if (col1 >= 0) {
            atom_1_comp_id = gemmi::cif::as_int(loop.val(j, col1));
        } else {
            ierr_tot++;
        }

        int col2 = loop.find_tag("atom_2_comp_id");
        if (col2 >= 0) {
            atom_2_comp_id = gemmi::cif::as_int(loop.val(j, col2));
        } else {
            ierr_tot++;
        }

        int col3 = loop.find_tag("atom_3_comp_id");
        if (col3 >= 0) {
            atom_3_comp_id = gemmi::cif::as_int(loop.val(j, col3));
        } else {
            ierr_tot++;
        }

        int col4 = loop.find_tag("value_angle");
        if (col4 >= 0) {
            value_angle = std::stod(loop.val(j, col4));
        } else {
            ierr_tot++;
        }

        int col5 = loop.find_tag("value_angle_esd");
        if (col5 >= 0) {
            value_angle_esd = std::stod(loop.val(j, col5));
        } else {
            ierr_tot++;
        }

        if (ierr_tot == 0) {
            link_add_angle_gemmi(link_id,
                atom_1_comp_id, atom_2_comp_id, atom_3_comp_id,
                atom_id_1, atom_id_2, atom_id_3,
                value_angle, value_angle_esd,
                dict_link_res_restraints);
        } else {
            std::cout << "problem reading link angle mmCIFLoop" << std::endl;
        }
    }
}

} // namespace protein_geometry
} // namespace coot