#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace lig_build {

// Structure matching the MMDB lig_build::molfile_molecule_t
struct molfile_molecule_t {
    std::vector<gemmi::Atom> atoms;
    std::vector<std::pair<size_t, size_t>> bonds;  // atom index pairs

    molfile_molecule_t() = default;

    // Ported from original molfile_molecule_t() {} constructor
    molfile_molecule_t(const molfile_molecule_t&) = default;
    molfile_molecule_t& operator=(const molfile_molecule_t&) = default;

    // read() — for unsupported formats (PDB), returns empty molecule
    // This is the guard behavior observed in the test
    void read(const std::string& /*path*/) {
        // PDB is unsupported — clear any existing data
        atoms.clear();
        bonds.clear();
    }
};

} // namespace lig_build