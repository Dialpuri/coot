#pragma once

#include <gemmi/cif.hpp>
#include <string>
#include <vector>
#include "coot/geometry/protein-geometry.hh"

namespace coot {

//! gemmi-port of protein_geometry::pdbx_chem_comp_description_generator
//! Takes a coot::protein_geometry reference and a gemmi cif::Loop,
//! populating description_generation entries in dict_res_restraints.
inline void pdbx_chem_comp_description_generator_gemmi(
    coot::protein_geometry& geom,
    gemmi::cif::Loop& loop,
    int imol_enc)
{
    int comp_id_col       = loop.find_tag("comp_id");
    int program_name_col  = loop.find_tag("program_name");
    int program_version_col = loop.find_tag("program_version");
    int descriptor_col    = loop.find_tag("descriptor");

    for (size_t j = 0; j < loop.length(); j++) {
        std::string comp_id;
        std::string program_name;
        std::string program_version;
        std::string descriptor;

        if (comp_id_col >= 0)         comp_id         = loop.val(j, comp_id_col);
        if (program_name_col >= 0)    program_name    = loop.val(j, program_name_col);
        if (program_version_col >= 0) program_version = loop.val(j, program_version_col);
        if (descriptor_col >= 0)      descriptor      = loop.val(j, descriptor_col);

        // Only proceed if comp_id tag exists (mirrors original ierr==0 guard)
        if (comp_id_col >= 0) {
            coot::pdbx_chem_comp_description_generator_t dg(
                program_name, program_version, descriptor);
            int idx = geom.get_monomer_restraints_index(comp_id, imol_enc, true);
            if (idx >= 0) {
                // operator[] is const-only; use const_cast to modify
                auto& entry = const_cast<
                    std::pair<int, dictionary_residue_restraints_t>&
                >(geom[idx]);
                entry.second.description_generation = dg;
            }
        }
    }
}

} // namespace coot