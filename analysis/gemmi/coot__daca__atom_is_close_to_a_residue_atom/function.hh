#pragma once
#include <gemmi/model.hpp>
#include <string>

namespace coot {
namespace daca {

// Check if an atom is close to any atom in a residue (within 4.9 Angstroms)
bool atom_is_close_to_a_residue_atom_gemmi(const gemmi::Atom& at, const gemmi::Residue& residue);

}
}