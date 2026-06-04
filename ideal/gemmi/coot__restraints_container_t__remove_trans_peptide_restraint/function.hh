#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <cstddef>

// Restrained type constant — matches the coot enum value from
// coot/ideal/refinement-results-mini-stats.hh
// (This is a standalone port; the constant is reproduced here to avoid
// pulling in the full coot headers.)
constexpr int TRANS_PEPTIDE_RESTRAINT_VAL = 2048;

struct simple_restraint_gemmi {
    int restraint_type;
    int atom_index_1;
    int atom_index_2;
    int atom_index_3;
    int atom_index_4;
    double value;
    double sigma;
    double weight;
    bool is_active;

    simple_restraint_gemmi()
        : restraint_type(0), atom_index_1(0), atom_index_2(0),
          atom_index_3(0), atom_index_4(0), value(0), sigma(1),
          weight(1), is_active(true) {}

    void close() { is_active = false; }
};

// gemmi port of restraints_container_t::remove_trans_peptide_restraint
//
// Scans all restraints for TRANS_PEPTIDE_RESTRAINT entries
// whose four atoms belong to the two given residues in the order
// (first, first, second, second), then calls close() on each matching
// restraint.
//
// Null arguments are handled gracefully (no-op) to match the original
// MMDB behavior.
//
// `atom_residues` is a parallel vector: atom_residues[i] is the
// gemmi::Residue* that owns the atom at global index i.
void remove_trans_peptide_restraint_gemmi(
    std::vector<simple_restraint_gemmi>& restraints_vec,
    const std::vector<gemmi::Residue*>& atom_residues,
    gemmi::Residue* first,
    gemmi::Residue* second);