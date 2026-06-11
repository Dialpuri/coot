#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <vector>
#include <string>

#include "coot/geometry/protein-geometry.hh"
#include "coot-utils/coot-coord-utils.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__protein_geometry__get_group__10c420/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_alt_confs/gemmi/function.hh"

namespace coot {

struct torsion_atom_quad_gemmi {
    gemmi::Atom *atom_1;
    gemmi::Atom *atom_2;
    gemmi::Atom *atom_3;
    gemmi::Atom *atom_4;
    double angle;
    double angle_esd;
    int period;
    std::string name;
    std::string residue_name;

    torsion_atom_quad_gemmi(gemmi::Atom *a1, gemmi::Atom *a2, gemmi::Atom *a3, gemmi::Atom *a4,
                            double ang, double esd, int per)
        : atom_1(a1), atom_2(a2), atom_3(a3), atom_4(a4),
          angle(ang), angle_esd(esd), period(per) {}
};

std::vector<coot::torsion_atom_quad_gemmi>
torsionable_bonds_monomer_internal_quads_gemmi(
    const gemmi::Residue &residue,
    std::vector<gemmi::Atom*> atom_selection,
    bool include_pyranose_ring_torsions_flag,
    coot::protein_geometry *geom_p);

} // namespace coot
