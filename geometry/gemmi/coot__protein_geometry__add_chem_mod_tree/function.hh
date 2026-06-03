#pragma once

#include "geometry/protein-geometry.hh"

// Lightweight mmCIF-loop-like structure for _chem_mod_tree rows.
// Mirrors the original mmdb::mmcif::PLoop interface just enough for the port.
struct ChemModTreeLoop {
    std::vector<std::string> tags;              // column names
    std::vector<std::vector<std::string>> rows; // one vector per row

    int GetLoopLength() const { return static_cast<int>(rows.size()); }

    const char* GetString(const std::string& tag, int row, int& ierr) const {
        ierr = 0;
        int col = -1;
        for (int i = 0; i < static_cast<int>(tags.size()); ++i) {
            if (tags[i] == tag) { col = i; break; }
        }
        if (col < 0 || row < 0 || row >= static_cast<int>(rows.size())) {
            ierr = 1; return nullptr;
        }
        if (rows[row].size() <= static_cast<std::size_t>(col)) {
            ierr = 1; return nullptr;
        }
        return rows[row][col].c_str();
    }
};

namespace coot {

// Gemmi port of protein_geometry::add_chem_mod_tree:
// reads chem_mod_tree loop data and populates the mods map.
void add_chem_mod_tree_gemmi(const ChemModTreeLoop& loop,
                             protein_geometry& geom);

} // namespace coot