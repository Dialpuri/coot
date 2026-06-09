#include "function.hh"

#include <gemmi/neighbor.hpp>
#include <gemmi/polyheur.hpp>

#include <iostream>
#include <chrono>
#include <set>

#include "geometry/protein-geometry.hh"
#include "geometry/residue-and-atom-specs.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

std::vector<h_bond_result>
make_h_bond_restraints_from_res_vec_auto_gemmi(
    const std::vector<gemmi::SeqId>& seq_ids,
    const std::vector<std::string>& res_names,
    const std::vector<std::string>& chain_ids,
    gemmi::Structure* st,
    const coot::protein_geometry& geom)
{
  auto tp_0 = std::chrono::high_resolution_clock::now();

  // Build atom index map: gemmi::Atom* -> int
  std::map<gemmi::Atom*, int> atom_index_map;

  // Collect all atoms from the given residues
  auto tp_1 = std::chrono::high_resolution_clock::now();
  for (unsigned int i = 0; i < seq_ids.size(); i++) {
    // residue_spec_t(chain_id, res_no, resname)
    residue_spec_t spec(chain_ids[i], static_cast<int>(seq_ids[i].num.value), res_names[i]);
    
    std::vector<gemmi::Atom*> atoms = coot::select_atoms_gemmi(spec, st);
    for (auto* a : atoms) {
      if (a) {
        atom_index_map[a] = static_cast<int>(atom_index_map.size());
      }
    }
  }
  auto tp_2 = std::chrono::high_resolution_clock::now();

  // Get all atoms that were selected
  std::vector<gemmi::Atom*> all_atoms;
  all_atoms.reserve(atom_index_map.size());
  for (auto& kv : atom_index_map) {
    all_atoms.push_back(kv.first);
  }

  // Build a set for O(n) lookup
  std::set<gemmi::Atom*> atom_set(all_atoms.begin(), all_atoms.end());

  // Use gemmi::NeighborSearch to find potential H-bonds
  gemmi::Model& model = st->models[0];
  gemmi::NeighborSearch ns(model, st->cell, 4.0); // max radius 4.0A for H-bonds
  ns.populate(true); // include hydrogens

  auto tp_3 = std::chrono::high_resolution_clock::now();

  // Track pairs we've already found (to avoid duplicates)
  // Use (min_index, max_index) as key
  std::set<std::pair<int,int>> seen_pairs;
  
  std::vector<h_bond_result> results;

  // For each atom, find neighbors that could form H-bonds
  for (auto* atom : all_atoms) {
    if (atom == nullptr) continue;
    
    // Find nearby atoms within 3.5 A (typical H-bond distance)
    std::vector<gemmi::NeighborSearch::Mark*> neighbors =
        ns.find_atoms(atom->pos, '\0', 0.0, 3.5);
    
    for (auto* mark : neighbors) {
      gemmi::CRA cra = mark->to_cra(model);
      gemmi::Atom* neighbor_atom = cra.atom;
      if (!neighbor_atom) continue;
      if (neighbor_atom == atom) continue;
      // Only consider atoms in our selection
      if (atom_set.find(neighbor_atom) == atom_set.end()) continue;

      // Use co_gemmi to get coordinates and compute distance
      clipper::Coord_orth pos_a = co_gemmi(atom);
      clipper::Coord_orth pos_b = co_gemmi(neighbor_atom);
      clipper::Coord_orth diff = pos_a - pos_b;
      double dist = sqrt(diff.lengthsq());
      
      if (dist < 1.0) continue; // too close, likely bonded
      
      // Check if this is a plausible H-bond: one atom should be N or O (donor/acceptor)
      // and the other should be N, O, or H
      std::string elem_a = atom->element.name();
      std::string elem_b = neighbor_atom->element.name();
      
      // Typical H-bond: donor (N or O with H attached) - acceptor (N or O with lone pair)
      // Simplified check: at least one should be N or O
      bool is_hbond_candidate = 
          (elem_a == "N" || elem_a == "O") && 
          (elem_b == "N" || elem_b == "O" || elem_b == "H");
      
      if (is_hbond_candidate) {
        // Avoid duplicates: only record (min_index, max_index) once
        int idx_a = atom_index_map[atom];
        int idx_b = atom_index_map[neighbor_atom];
        int min_idx = std::min(idx_a, idx_b);
        int max_idx = std::max(idx_a, idx_b);
        
        auto key = std::make_pair(min_idx, max_idx);
        if (seen_pairs.find(key) == seen_pairs.end()) {
          seen_pairs.insert(key);
          
          h_bond_result hb;
          hb.donor = atom;
          hb.acceptor = neighbor_atom;
          hb.dist = dist;
          results.push_back(hb);
        }
      }
    }
  }

  auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
  auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
  auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
  std::cout << "------------------- timing: " << d10 << " " << d21 << " " << d32
            << " milliseconds to find " << results.size() << " H-bonds " << std::endl;
  std::cout << "DEBUG:: made " << results.size() << " hydrogen bonds " << std::endl;

  return results;
}

} // namespace coot
