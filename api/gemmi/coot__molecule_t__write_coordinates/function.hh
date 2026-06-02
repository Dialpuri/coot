#pragma once

#include <gemmi/model.hpp>
#include <gemmi/to_mmcif.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/to_pdb.hpp>
#include <fstream>
#include <string>
#include <sstream>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__correct_link_distances/gemmi/function.hh"

namespace coot {

// ── helper utilities (no gemmi ports exist) ──────────────────────

inline std::string file_name_extension(const std::string& file_name) {
    auto pos = file_name.rfind('.');
    if (pos == std::string::npos) return "";
    return file_name.substr(pos);
}

inline bool extension_is_for_shelx_coords(const std::string& ext) {
    return (ext == ".ins" || ext == ".res" || ext == ".xyz");
}

// ── inline write_coords_pdb (gemmi version, mirrors generated port) ─

inline int write_coords_pdb_gemmi(gemmi::Structure& st, const std::string& file_name) {
    // Mirrors coot::write_coords_pdb_gemmi from generated port:
    //   1. correct link distances (cis-peptide cleanup is no-op for gemmi)
    //   2. write PDB
    coot::util::correct_link_distances_gemmi(st);

    std::ofstream ofs(file_name);
    if (!ofs) return -1;
    gemmi::write_pdb(st, ofs);
    return 0;
}

// ── port of coot::molecule_t::write_coordinates ──────────────────

// Original:
//   int coot::molecule_t::write_coordinates(const std::string &file_name) const
//
// gemmi port: takes a gemmi::Structure + file_name + n_selected_atoms
inline int write_coordinates_gemmi(
    gemmi::Structure& st,
    const std::string& file_name,
    int n_selected_atoms
) {
    int err = 1;
    if (n_selected_atoms > 0) {
        std::string ext = file_name_extension(file_name);
        if (extension_is_for_shelx_coords(ext)) {
            // shelx — no port for write_shelx_ins_file; leave err=1
        } else if (ext == ".cif") {
            gemmi::cif::Document doc = gemmi::make_mmcif_document(st);
            std::ofstream ofs(file_name);
            if (!ofs) {
                err = -1;
            } else {
                gemmi::cif::write_cif_to_stream(ofs, doc);
                err = 0;
            }
        } else {
            err = write_coords_pdb_gemmi(st, file_name);
        }
    }
    return err;
}

} // namespace coot