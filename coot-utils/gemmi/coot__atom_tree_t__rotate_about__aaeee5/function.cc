#include "function.hh"

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <vector>
#include <string>
#include <utility>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "coot-utils/map-index.hh"
#include "coot-utils/atom-vertex.hh"
#include "utils/coot-utils.hh"
#include "clipper/core/clipper_util.h"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_tree_t__rotate_internal/gemmi/function.hh"

namespace coot {

/* ---- helpers (no _gemmi ports exist yet) ---- */

static std::vector<map_index_t> uniquify_atom_indices(
    const std::vector<map_index_t>& indices) {
  std::vector<map_index_t> result;
  for (const auto& idx : indices) {
    bool found = false;
    for (const auto& e : result) {
      if (e.index() == idx.index()) { found = true; break; }
    }
    if (!found) result.push_back(idx);
  }
  return result;
}

static std::pair<int, std::vector<map_index_t>> get_forward_atoms(
    const std::vector<atom_vertex>& atom_vertex_vec,
    const std::vector<map_index_t>& atom_indices,
    map_index_t start_index,
    map_index_t end_index) {
  int n_called = 0;
  std::vector<map_index_t> fwd;
  std::vector<int> stack;

  int si = start_index.index();
  if (si >= 0 && static_cast<size_t>(si) < atom_vertex_vec.size()) {
    for (int f : atom_vertex_vec[si].forward) {
      if (f != end_index.index()) stack.push_back(f);
    }
  }

  while (!stack.empty()) {
    int cur = stack.back(); stack.pop_back();
    n_called++;
    if (cur >= 0 && static_cast<size_t>(cur) < atom_indices.size()) {
      fwd.push_back(atom_indices[cur]);
    }
    if (cur >= 0 && static_cast<size_t>(cur) < atom_vertex_vec.size()) {
      for (int f : atom_vertex_vec[cur].forward) {
        if (f != end_index.index()) stack.push_back(f);
      }
    }
  }
  return {n_called, fwd};
}

static std::vector<map_index_t> complementary_indices(
    const std::vector<map_index_t>& all_atom_indices,
    map_index_t idx2,
    map_index_t idx3,
    const std::vector<map_index_t>& moving) {
  std::vector<map_index_t> result;
  for (const auto& a : all_atom_indices) {
    int i = a.index();
    if (i == idx2.index() || i == idx3.index()) continue;
    bool in_moving = false;
    for (const auto& m : moving) {
      if (m.index() == i) { in_moving = true; break; }
    }
    if (!in_moving) result.push_back(a);
  }
  return result;
}

static double quad_to_torsion(
    const std::vector<atom_vertex>& atom_vertex_vec,
    const std::vector<const gemmi::Atom*>& atom_selection,
    int index2) {
  if (static_cast<size_t>(index2) >= atom_vertex_vec.size()) return 0.0;
  const auto& tq = atom_vertex_vec[index2].torsion_quad;
  if (!tq.first) return 0.0;
  const auto& quad = tq.second;
  if (atom_selection.empty()) return 0.0;

  auto get_atom = [&](int qi) -> const gemmi::Atom* {
    if (qi < 0 || qi >= static_cast<int>(atom_selection.size())) return nullptr;
    return atom_selection[qi];
  };

  // atom_index_quad has fields: index1, index2, index3, index4
  const gemmi::Atom* a0 = get_atom(quad.index1);
  const gemmi::Atom* a1 = get_atom(quad.index2);
  const gemmi::Atom* a2 = get_atom(quad.index3);
  const gemmi::Atom* a3 = get_atom(quad.index4);
  if (!a0 || !a1 || !a2 || !a3) return 0.0;

  gemmi::Vec3 p0(a0->pos), p1(a1->pos), p2(a2->pos), p3(a3->pos);
  gemmi::Vec3 v01 = p1 - p0;
  gemmi::Vec3 v12 = p2 - p1;
  gemmi::Vec3 v23 = p3 - p2;

  gemmi::Vec3 n1 = v01.cross(v12);
  gemmi::Vec3 n2 = v12.cross(v23);
  gemmi::Vec3 p  = n1.cross(v12);

  double l1 = n1.length();
  double l2 = n2.length();
  double lp = p.length();
  if (l1 < 1e-10 || l2 < 1e-10 || lp < 1e-10) return 0.0;

  double cos_val = n1.dot(n2) / (l1 * l2);
  double sin_val = p.dot(v12) / (l1 * l2 * v12.length());
  return std::atan2(sin_val, cos_val);
}

/* ---- main function ---- */

double atom_tree_t_rotate_about_gemmi(
    const std::vector<atom_vertex>& atom_vertex_vec,
    const std::vector<const gemmi::Atom*>& atom_selection,
    gemmi::Residue* residue,
    int index2,
    int index3,
    double angle,
    bool reversed_flag) {

  bool debug = false;
  double new_torsion = 0.0;

  // Validate indices
  if (index2 == -1) {
    throw std::runtime_error("Ooops! rotate_about() Bad atom index: index-2");
  }
  if (index3 == -1) {
    throw std::runtime_error("Ooops! rotate_about() Bad atom index: index-3");
  }

  bool index2_has_fwd = (index2 >= 0 && static_cast<size_t>(index2) < atom_vertex_vec.size()
                         && atom_vertex_vec[index2].forward.size() > 0);
  bool index3_has_fwd = (index3 >= 0 && static_cast<size_t>(index3) < atom_vertex_vec.size()
                         && atom_vertex_vec[index3].forward.size() > 0);

  if (!index2_has_fwd && !index3_has_fwd) {
    std::string s = "Neither index2 ";
    s += coot::util::int_to_string(index2);
    s += " nor index3 ";
    s += coot::util::int_to_string(index3);
    s += " has forward atoms!";
    throw std::runtime_error(s);
  }

  bool internal_reversed = false;
  bool index3_is_forward = false;

  // Check if index3 is in the forward atoms of index2
  if (index2 >= 0 && static_cast<size_t>(index2) < atom_vertex_vec.size()) {
    for (unsigned int ifo = 0; ifo < atom_vertex_vec[index2].forward.size(); ifo++) {
      if (atom_vertex_vec[index2].forward[ifo] == index3) {
        index3_is_forward = true;
        break;
      }
    }
  }

  if (!index3_is_forward) {
    bool index2_is_forward = false;
    if (index3 >= 0 && static_cast<size_t>(index3) < atom_vertex_vec.size()) {
      for (unsigned int ifo = 0; ifo < atom_vertex_vec[index3].forward.size(); ifo++) {
        if (atom_vertex_vec[index3].forward[ifo] == index2) {
          index2_is_forward = true;
          break;
        }
      }
    }

    if (index2_is_forward) {
      if (debug)
        std::cout << "    index2 was the forward atom of index3 - swapping and setting internal_reversed" << std::endl;
      std::swap(index2, index3);
      index3_is_forward = true;
      internal_reversed = true;
    }
  }

  if (index3_is_forward) {
    // Build atom_indices from atom_selection
    std::vector<map_index_t> atom_indices;
    for (size_t i = 0; i < atom_selection.size(); i++) {
      if (atom_selection[i]) atom_indices.push_back(map_index_t(static_cast<int>(i)));
    }

    auto [n_calls, moving_atom_indices] = get_forward_atoms(
        atom_vertex_vec, atom_indices,
        map_index_t(index3), map_index_t(index3));

    std::vector<map_index_t> unique_moving_atom_indices =
        uniquify_atom_indices(moving_atom_indices);

    bool xor_reverse = reversed_flag ^ internal_reversed;
    if (xor_reverse) {
      unique_moving_atom_indices = complementary_indices(
          atom_indices, map_index_t(index2), map_index_t(index3),
          unique_moving_atom_indices);
    }

    const gemmi::Atom* at2 = nullptr;
    const gemmi::Atom* at3 = nullptr;

    if (!atom_selection.empty()) {
      if (index2 >= 0 && index2 < static_cast<int>(atom_selection.size()))
        at2 = atom_selection[index2];
      if (index3 >= 0 && index3 < static_cast<int>(atom_selection.size()))
        at3 = atom_selection[index3];
    }
    if (residue) {
      int n_res = static_cast<int>(residue->atoms.size());
      if (index2 >= 0 && index2 < n_res) at2 = &residue->atoms[index2];
      if (index3 >= 0 && index3 < n_res) at3 = &residue->atoms[index3];
    }

    if (at2 && at3) {
      gemmi::Vec3 base_atom_pos(at2->pos);
      gemmi::Vec3 third_atom(at3->pos);
      gemmi::Vec3 direction = third_atom - base_atom_pos;
      if (xor_reverse) {
        direction = base_atom_pos - third_atom;
        base_atom_pos = third_atom;
      }

      // Create a mutable copy of atoms for rotate_internal_gemmi
      std::vector<gemmi::Atom> atoms_copy;
      if (residue) {
        for (const auto& a : residue->atoms) {
          atoms_copy.push_back(a);
        }
      }

      // Call the existing gemmi port of rotate_internal
      coot::atom_tree_t dummy_tree;
      dummy_tree.rotate_internal_gemmi(
          atoms_copy,
          unique_moving_atom_indices, direction, base_atom_pos,
          clipper::Util::d2rad(angle));

      // Copy positions back to the actual atoms in the structure
      if (residue && atoms_copy.size() == residue->atoms.size()) {
        for (size_t i = 0; i < residue->atoms.size(); i++) {
          residue->atoms[i].pos = atoms_copy[i].pos;
        }
      }

      // Compute new_torsion if torsion quad is set
      if (static_cast<size_t>(index2) < atom_vertex_vec.size() &&
          atom_vertex_vec[index2].torsion_quad.first) {
        new_torsion = quad_to_torsion(atom_vertex_vec, atom_selection, index2);
      }
    } else {
      std::cout << "ERROR:: null atom rotate_about() - this should not happen" << std::endl;
    }
  }
  return new_torsion;
}

} // namespace coot