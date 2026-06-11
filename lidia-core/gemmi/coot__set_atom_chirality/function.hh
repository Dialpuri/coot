#pragma once

// Re-export the verified 4-arg gemmi port of coot::set_atom_chirality.
// This header pulls in the 4-arg overload (const gemmi::Atom*, const gemmi::Residue*)
// and transitively includes the 2-arg overload (const coot::dict_atom&) as well.
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__set_atom_chirality__4e3f81/gemmi/function.hh"
