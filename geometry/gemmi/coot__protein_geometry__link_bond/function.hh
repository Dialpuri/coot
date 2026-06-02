#pragma once

#include <gemmi/cifdoc.hpp>
#include <string>
#include <vector>
#include <iostream>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__protein_geometry__link_add_bond/gemmi/function.hh"

namespace coot {
namespace protein_geometry {

int link_bond_gemmi(gemmi::cif::Loop& loop, std::vector<dictionary_residue_link_restraints_t>& dict_link_res_restraints) {
    std::string link_id;
    std::string atom_id_1, atom_id_2;
    double value_dist, value_dist_esd;
    int atom_1_comp_id = 0, atom_2_comp_id = 0;

    int n_link_bonds = 0;

    for (size_t j = 0; j < loop.length(); j++) {
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

        int col3 = loop.find_tag("value_dist");
        if (col3 >= 0) {
            value_dist = std::stod(loop.val(j, col3));
        } else {
            ierr_tot++;
        }

        int col4 = loop.find_tag("value_dist_esd");
        if (col4 >= 0) {
            value_dist_esd = std::stod(loop.val(j, col4));
        } else {
            ierr_tot++;
        }

        if (ierr_tot == 0) {
            coot::link_add_bond_gemmi(link_id,
                atom_1_comp_id, atom_2_comp_id,
                atom_id_1, atom_id_2,
                value_dist, value_dist_esd,
                dict_link_res_restraints);
            n_link_bonds++;
        } else {
            std::cout << "problem reading bond mmCIFLoop" << std::endl;
        }
    }

    return n_link_bonds;
}

} // namespace protein_geometry
} // namespace coot