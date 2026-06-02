#pragma once
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <gemmi/model.hpp>
#include <clipper/core/clipper_util.h>
#include "ccp4mg-utils/cartesian.h"
#include "ccp4mg-utils/mgtree.h"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__chi_angles__get_atom_index_pairs/gemmi/function.hh"
#include "geometry/residue-and-atom-specs.hh"

namespace coot {
namespace chi_angles {

std::pair<short int, float>
change_by_internal_gemmi(int ichi,
                         double diff,
                         const std::vector<coot::chi_angles::atom_name_pair>& atom_name_pairs,
                         const std::vector<std::vector<int>>& contact_indices,
                         gemmi::Residue& residue,
                         const coot::atom_spec_t& tree_base_atom) {

    std::pair<short int, float> p(0, 0.0);

    int nResidueAtoms = static_cast<int>(residue.atoms.size());

    // Build coordinates vector
    std::vector<::Cartesian> coords;
    for (int i = 0; i < nResidueAtoms; i++) {
        ::Cartesian c(residue.atoms[i].pos.x,
                      residue.atoms[i].pos.y,
                      residue.atoms[i].pos.z);
        coords.push_back(c);
    }

    // Get atom index pairs using the gemmi port
    std::vector<coot::chi_angles::atom_index_pair> atom_index_pairs =
        coot::chi_angles::get_atom_index_pairs_gemmi(atom_name_pairs, residue);

    // Find base index
    int base_index = 0;
    if (tree_base_atom.int_user_data != -999) {
        for (int i = 0; i < nResidueAtoms; i++) {
            if (tree_base_atom.atom_name == residue.atoms[i].name) {
                base_index = i;
                break;
            }
        }
    }

    Tree tree;
    tree.SetCoords(coords, base_index, contact_indices);

    int ibond = ichi - 1;
    if (ibond >= static_cast<int>(atom_index_pairs.size())) {
        std::cout << "ERROR: impossible ibond " << ibond
                  << " (wanted rotamer index " << ichi << ") in "
                  << atom_index_pairs.size() << " rotamer atom pairs" << std::endl;
        p.first = 2;
        return p;
    }

    float tors_orig;
    float tors;
    TreeVertex* tv = tree.GetCoord(atom_index_pairs[ibond].index2);
    if (tv->GetNumberOfChildren() > 0) {
        TreeVertex* tvc0 = tv->GetChild(0);
        tors_orig = tvc0->GetParentDihedralAngle();

        tors = clipper::Util::d2rad(diff);
        p.second = tors_orig + tors; // in radians

        tree.RotateAboutBond(atom_index_pairs[ibond].index2,
                             atom_index_pairs[ibond].index1, tors);

        std::vector<::Cartesian> coords_rotated = tree.GetAllCartesians();

        if (static_cast<int>(coords_rotated.size()) != nResidueAtoms) {
            std::cout << "disaster in atom selection, trees, dunbrack\n";
            p.first = 3;
        } else {
            for (int iat = 0; iat < nResidueAtoms; iat++) {
                residue.atoms[iat].pos.x = coords_rotated[iat].get_x();
                residue.atoms[iat].pos.y = coords_rotated[iat].get_y();
                residue.atoms[iat].pos.z = coords_rotated[iat].get_z();
            }
        }
    } else {
        std::cout << "WARNING: this vertex " << atom_index_pairs[ibond].index2
                  << " has no children (strangely)\n";
        std::cout << "         bond: " << ibond << " index2: "
                  << atom_index_pairs[ibond].index2 << "\n";
        TreeVertex* tv1 = tree.GetCoord(atom_index_pairs[ibond].index1);
        std::cout << "         tv1 (index: " << atom_index_pairs[ibond].index1
                  << ") has " << tv1->GetNumberOfChildren()
                  << " children \n";
        p.first = 1;
    }

    return p;
}

} // namespace chi_angles
} // namespace coot