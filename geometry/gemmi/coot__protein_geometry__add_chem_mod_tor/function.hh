#pragma once

#include <string>
#include <gemmi/cif.hpp>
#include "geometry/protein-geometry.hh"

void add_chem_mod_tor_gemmi(
    coot::protein_geometry& geom,
    const gemmi::cif::Loop& cif_loop);