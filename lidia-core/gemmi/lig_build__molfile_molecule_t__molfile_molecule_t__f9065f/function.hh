#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>

namespace lig_build {

// Gemmi port of molfile_atom_t — stores atom position, element, and name
struct molfile_atom_t {
    gemmi::Position pos;
    std::string element;
    std::string atom_name;

    molfile_atom_t(gemmi::Position p, std::string elem, std::string name)
        : pos(std::move(p)), element(std::move(elem)), atom_name(std::move(name)) {}
};

// Gemmi port of molfile_molecule_t — holds atoms and bonds
struct molfile_molecule_t {
    std::vector<gemmi::Atom> atoms;
    std::vector<std::pair<size_t, size_t>> bonds;

    // Default constructor (ported from MMDB molfile_molecule_t() {})
    molfile_molecule_t() = default;

    // Ported constructor — factory returning an empty molecule
    static molfile_molecule_t molfile_molecule_t_gemmi() {
        return molfile_molecule_t{};
    }

    // Add an atom to the molecule
    void add_molfile_atom(const molfile_atom_t& atom);
};

inline void molfile_molecule_t::add_molfile_atom(const molfile_atom_t& atom) {
    gemmi::Atom a;
    a.pos     = atom.pos;
    a.element = gemmi::Element(atom.element);
    a.name    = atom.atom_name;
    atoms.push_back(a);
}

} // namespace lig_build