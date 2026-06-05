#pragma once

#include <gemmi/model.hpp>
#include <gemmi/symmetry.hpp>
#include <vector>
#include <string>
#include <set>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__water_coordination_t__add_contacts/gemmi/function.hh"

namespace coot {
namespace util {

inline void init_internal_gemmi(
    water_coordination_t& wc,
    gemmi::Structure* st,
    double radius,
    bool do_metals_only_flag)
{
  if (!st) return;

  gemmi::Model* mol = &st->first_model();

  double min_dist = 0.5;
  double max_dist = radius;

  // Metal elements list from the original MMDB code
  static const std::set<std::string> metals = {
    "MG","CA","K","NA","LI","RB","BE","BA","FR","CS",
    "SR","RA","SC","TI","V","CR","MN","FE","CO","NI",
    "CU","ZN","ZR","NB","MO","RU","RH","AG","CD","W",
    "OS","IR","PT","AU","HG"
  };

  // Build water selection: HOH residues (unless metals-only)
  std::vector<gemmi::Atom const*> water_selection;
  if (!do_metals_only_flag) {
    for (gemmi::Chain& chain : mol->chains) {
      for (gemmi::Residue& res : chain.residues) {
        if (res.name == "HOH") {
          for (gemmi::Atom& atom : res.atoms) {
            water_selection.push_back(&atom);
          }
        }
      }
    }
  }

  // Always add metal atoms to water selection
  for (gemmi::Chain& chain : mol->chains) {
    for (gemmi::Residue& res : chain.residues) {
      for (gemmi::Atom& atom : res.atoms) {
        if (metals.count(atom.element.name())) {
          water_selection.push_back(&atom);
        }
      }
    }
  }

  // Build all-atoms selection (everything in the model)
  std::vector<gemmi::Atom const*> all_selection;
  for (gemmi::Chain& chain : mol->chains) {
    for (gemmi::Residue& res : chain.residues) {
      for (gemmi::Atom& atom : res.atoms) {
        all_selection.push_back(&atom);
      }
    }
  }

  // Identity 4x4 transformation matrix
  double my_matt[4][4] = {
    {1.0, 0.0, 0.0, 0.0},
    {0.0, 1.0, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0},
    {0.0, 0.0, 0.0, 1.0}
  };

  // Check for non-trivial symmetry
  const gemmi::SpaceGroup* sg = st->find_spacegroup();
  gemmi::GroupOps group_ops = sg ? sg->operations() : gemmi::GroupOps{};
  bool has_symmetry = (sg != nullptr && sg->number != 1 && !group_ops.sym_ops.empty());

  if (!has_symmetry) {
    // No symmetry: single call with identity matrix
    wc.add_contacts_gemmi(water_selection, all_selection, min_dist, max_dist, my_matt);
  } else {
    // Has symmetry: iterate over space group ops and lattice translations
    for (int ix = -1; ix <= 1; ix++) {
      for (int iy = -1; iy <= 1; iy++) {
        for (int iz = -1; iz <= 1; iz++) {
          for (const gemmi::Op& op : group_ops.sym_ops) {
            // Build 4x4 matrix from gemmi::Op using float_seitz()
            auto seitz = op.float_seitz();
            for (int i = 0; i < 4; i++)
              for (int j = 0; j < 4; j++)
                my_matt[i][j] = seitz[i][j];
            wc.add_contacts_gemmi(water_selection, all_selection, min_dist, max_dist, my_matt);
          }
        }
      }
    }
  }
}

} // namespace util
} // namespace coot