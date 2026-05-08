#ifndef COOT_MOLECULE_APPEND_TO_MOLECULE_GEMMI_HH
#define COOT_MOLECULE_APPEND_TO_MOLECULE_GEMMI_HH

#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

// Port of coot::molecule_t::append_to_molecule to gemmi
// Appends atoms from water_mol to dest structure
// Returns 0 if new atoms added, 1 if chain already existed
int append_to_molecule_gemmi(gemmi::Structure& dest, const gemmi::Structure& water_mol);

}

#endif // COOT_MOLECULE_APPEND_TO_MOLECULE_GEMMI_HH