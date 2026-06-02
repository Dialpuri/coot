#pragma once

#include <gemmi/model.hpp>
#include <string>
#include "lidia-core/use-rdkit.hh"

namespace coot {
    class dictionary_residue_restraints_t;

    RDKit::RWMol rdkit_mol_gemmi(
        gemmi::CRA cra,
        const coot::dictionary_residue_restraints_t &restraints,
        const std::string &alt_conf,
        bool do_undelocalize);
}