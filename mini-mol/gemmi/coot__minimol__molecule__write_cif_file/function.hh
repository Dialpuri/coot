#pragma once

#include <gemmi/model.hpp>
#include <gemmi/to_mmcif.hpp>
#include <gemmi/to_cif.hpp>
#include <string>
#include <fstream>

namespace coot {
namespace minimol {
namespace molecule {

inline int write_cif_file_gemmi(
    const gemmi::Structure& st,
    const std::string& cif_filename)
{
    gemmi::cif::Document doc = gemmi::make_mmcif_document(st);
    std::ofstream out(cif_filename);
    if (!out) return 1;
    gemmi::cif::write_cif_to_stream(out, doc);
    out.close();
    return 0;
}

} // namespace molecule
} // namespace minimol
} // namespace coot