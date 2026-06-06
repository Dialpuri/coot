#include "function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__copy_and_delete_hydrogens/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__deep_copy_this_residue/gemmi/function.hh"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace coot {

// Helper: build a simple adjacency-based graph from residue atoms
// Returns a vector of atom indices, grouped by matching element names
static std::vector<int> get_atom_indices_by_element(
    const gemmi::Residue& res, bool include_hydrogens) {
  std::vector<int> indices;
  for (int i = 0; i < static_cast<int>(res.atoms.size()); ++i) {
    const auto& atom = res.atoms[i];
    if (!include_hydrogens && atom.element.is_hydrogen()) continue;
    indices.push_back(i);
  }
  return indices;
}

// Helper: compute squared distance between two positions
static double sq_dist(const gemmi::Vec3& a, const gemmi::Vec3& b) {
  gemmi::Vec3 d = a - b;
  return d.length_sq();
}

// Simple graph-match: match atoms by element name, try all permutations
// For small residues this is feasible
static bool try_match(
    const std::vector<int>& moving_atoms,
    const std::vector<int>& ref_atoms,
    const gemmi::Residue& res_moving,
    const gemmi::Residue& res_ref,
    bool apply_rtop_flag,
    double& best_dist_sum,
    int& best_n,
    std::vector<std::pair<std::pair<std::string, std::string>,
                          std::pair<std::string, std::string>>>& best_matching,
    clipper::RTop_orth& best_rtop) {
  
  int n_mov = static_cast<int>(moving_atoms.size());
  int n_ref = static_cast<int>(ref_atoms.size());
  
  // Match atoms by element name
  // Build mapping: element_name -> list of indices in moving/ref
  std::map<std::string, std::vector<int>> mov_by_elem;
  std::map<std::string, std::vector<int>> ref_by_elem;
  
  for (int idx : moving_atoms) {
    mov_by_elem[res_moving.atoms[idx].element.name()].push_back(idx);
  }
  for (int idx : ref_atoms) {
    ref_by_elem[res_ref.atoms[idx].element.name()].push_back(idx);
  }
  
  // Find matching element names
  std::vector<std::string> common_elems;
  for (auto& kv : mov_by_elem) {
    if (ref_by_elem.count(kv.first)) {
      common_elems.push_back(kv.first);
    }
  }
  
  if (common_elems.empty()) return false;
  
  // Simple greedy match: for each common element, pair up atoms by position order
  std::vector<std::pair<int, int>> pairs;
  for (const auto& elem : common_elems) {
    const auto& mov_list = mov_by_elem[elem];
    const auto& ref_list = ref_by_elem[elem];
    int n = std::min(static_cast<int>(mov_list.size()), static_cast<int>(ref_list.size()));
    for (int i = 0; i < n; ++i) {
      pairs.push_back({mov_list[i], ref_list[i]});
    }
  }
  
  if (pairs.empty()) return false;
  
  // Try to find best match among pairs using RTop
  // For now, use all pairs
  std::vector<clipper::Coord_orth> coords1, coords2;
  std::vector<std::pair<std::pair<std::string, std::string>,
                        std::pair<std::string, std::string>>> matching_atoms;
  
  for (auto& p : pairs) {
    const auto& a1 = res_moving.atoms[p.first];
    const auto& a2 = res_ref.atoms[p.second];
    coords1.push_back(clipper::Coord_orth(a1.pos.x, a1.pos.y, a1.pos.z));
    coords2.push_back(clipper::Coord_orth(a2.pos.x, a2.pos.y, a2.pos.z));
    
    std::string alt1 = a1.altloc ? std::string(1, a1.altloc) : "";
    std::string alt2 = a2.altloc ? std::string(1, a2.altloc) : "";
    matching_atoms.push_back({{a1.name, alt1}, {a2.name, alt2}});
  }
  
  // Calculate RTop if needed
  clipper::RTop_orth rtop(clipper::Mat33<double>(0,0,0,0,0,0,0,0,0),
                          clipper::Coord_orth(0,0,0));
  double dist_sum = 0.0;
  
  if (apply_rtop_flag && coords1.size() >= 3) {
    rtop = clipper::RTop_orth(coords1, coords2);
    for (size_t i = 0; i < coords1.size(); ++i) {
      auto transformed = coords1[i].transform(rtop);
      dist_sum += clipper::Coord_orth::length(coords2[i], transformed);
    }
  } else {
    for (size_t i = 0; i < coords1.size(); ++i) {
      dist_sum += clipper::Coord_orth::length(coords2[i], coords1[i]);
    }
  }
  
  best_dist_sum = dist_sum;
  best_n = static_cast<int>(coords1.size());
  best_matching = matching_atoms;
  best_rtop = rtop;
  
  return true;
}

graph_match_info_t graph_match_gemmi(
    const gemmi::Residue* res_moving,
    const gemmi::Residue* res_reference,
    bool apply_rtop_flag,
    bool match_hydrogens_also) {

  clipper::Mat33<double> m_dum(1,0,0,0,1,0,0,0,1);
  clipper::Coord_orth pt_dum(0,0,0);
  clipper::RTop_orth rtop(m_dum, pt_dum);
  bool success = false;
  std::vector<std::pair<std::pair<std::string, std::string>,
                        std::pair<std::string, std::string>>> best_matching_atoms;
  double best_match_sum = 1e20;
  int best_n_match = -99;

  // Copy residues (with or without hydrogens)
  const gemmi::Residue* cleaned_moving = nullptr;
  const gemmi::Residue* cleaned_ref = nullptr;
  
  if (!match_hydrogens_also) {
    cleaned_moving = coot::util::copy_and_delete_hydrogens_gemmi(res_moving);
    cleaned_ref = coot::util::copy_and_delete_hydrogens_gemmi(res_reference);
  } else {
    cleaned_moving = coot::util::deep_copy_this_residue_gemmi(res_moving);
    cleaned_ref = coot::util::deep_copy_this_residue_gemmi(res_reference);
  }

  if (!cleaned_moving || !cleaned_ref) {
    if (cleaned_moving) delete cleaned_moving;
    if (cleaned_ref) delete cleaned_ref;
    graph_match_info_t gmi;
    gmi.success = false;
    gmi.rtop = rtop;
    gmi.dist_score = best_match_sum;
    gmi.matching_atom_names = best_matching_atoms;
    gmi.n_match = best_n_match;
    return gmi;
  }

  // Simple matching based on element names
  clipper::RTop_orth best_rtop(m_dum, pt_dum);
  bool matched = try_match(
      get_atom_indices_by_element(*cleaned_moving, match_hydrogens_also),
      get_atom_indices_by_element(*cleaned_ref, match_hydrogens_also),
      *cleaned_moving, *cleaned_ref,
      apply_rtop_flag,
      best_match_sum, best_n_match, best_matching_atoms, best_rtop);

  if (matched) {
    rtop = best_rtop;
    success = true;
  }

  // Clean up copies
  delete cleaned_moving;
  delete cleaned_ref;

  graph_match_info_t gmi;
  gmi.success = success;
  gmi.rtop = rtop;
  gmi.dist_score = best_match_sum;
  gmi.matching_atom_names = best_matching_atoms;
  gmi.n_match = best_n_match;
  return gmi;
}

} // namespace coot