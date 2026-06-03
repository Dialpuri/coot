#pragma once

#include <coot/geometry/protein-geometry.hh>
#include <rdkit/GraphMol/RDKitBase.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/RDGeneral/Invariant.h>

#include <iostream>
#include <cstdio>
#include <utility>
#include <vector>
#include <string>

namespace coot {

void add_chem_comp_sp2_N_planes_gemmi(const RDKit::ROMol &mol, dictionary_residue_restraints_t *restraints);

}