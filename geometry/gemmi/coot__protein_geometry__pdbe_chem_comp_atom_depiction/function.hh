#pragma once

#include <gemmi/cif.hpp>
#include <string>
#include <vector>
#include <set>
#include "coot/geometry/protein-geometry.hh"

namespace coot {

//! gemmi-port of protein_geometry::pdbe_chem_comp_atom_depiction
//! Takes a coot::protein_geometry reference and a gemmi cif::Loop,
//! populating depiction entries in dict_res_restraints.
inline void pdbe_chem_comp_atom_depiction_gemmi(
    coot::protein_geometry& geom,
    gemmi::cif::Loop& loop,
    int imol_enc)
{
    int comp_id_col   = loop.find_tag("comp_id");
    int atom_id_col   = loop.find_tag("atom_id");
    int element_col   = loop.find_tag("element");
    int x_col         = loop.find_tag("model_Cartn_x");
    int y_col         = loop.find_tag("model_Cartn_y");
    int ordinal_col   = loop.find_tag("pdbx_ordinal");

    std::vector<coot::depiction_atom_t> dav;
    std::set<std::string> comp_id_set;

    for (size_t j = 0; j < loop.length(); j++) {
        std::string atom_id, element;
        double model_Cartn_x = 0.0, model_Cartn_y = 0.0;
        int pdbx_ordinal = 0;

        int ierr     = 0;
        int ierr_x   = 0;
        int ierr_y   = 0;
        int ierr_ord = 0;

        // comp_id
        if (comp_id_col >= 0) {
            std::string comp_id = loop.val(j, comp_id_col);
            comp_id_set.insert(comp_id);
        } else {
            ierr = 1;
        }

        // atom_id
        if (atom_id_col >= 0) {
            atom_id = loop.val(j, atom_id_col);
        } else {
            ierr = 1;
        }

        // element
        if (element_col >= 0) {
            element = loop.val(j, element_col);
        } else {
            ierr = 1;
        }

        // model_Cartn_x
        if (x_col >= 0) {
            try {
                model_Cartn_x = std::stod(loop.val(j, x_col));
            } catch (...) {
                ierr_x = 1;
            }
        } else {
            ierr_x = 1;
        }

        // model_Cartn_y
        if (y_col >= 0) {
            try {
                model_Cartn_y = std::stod(loop.val(j, y_col));
            } catch (...) {
                ierr_y = 1;
            }
        } else {
            ierr_y = 1;
        }

        // pdbx_ordinal
        if (ordinal_col >= 0) {
            try {
                pdbx_ordinal = std::stoi(loop.val(j, ordinal_col));
            } catch (...) {
                ierr_ord = 1;
            }
        } else {
            ierr_ord = 1;
        }

        if (ierr == 0 && ierr_x == 0 && ierr_y == 0 && ierr_ord == 0) {
            coot::depiction_atom_t da(atom_id, element, model_Cartn_x, model_Cartn_y, pdbx_ordinal);
            dav.push_back(da);
        }
    }

    if (!dav.empty()) {
        if (comp_id_set.size() == 1) {
            std::string comp_id = *comp_id_set.begin();
            coot::chem_comp_atom_depiction_t d(comp_id, dav);
            int idx = geom.get_monomer_restraints_index(comp_id, imol_enc, true);
            if (idx >= 0) {
                auto& entry = const_cast<
                    std::pair<int, dictionary_residue_restraints_t>&
                >(geom[idx]);
                entry.second.depiction = d;
            }
        }
    }
}

} // namespace coot