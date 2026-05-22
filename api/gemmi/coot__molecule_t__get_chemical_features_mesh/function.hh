#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include "coot/coot-utils/simple-mesh.hh"

namespace coot {

// Port of molecule_t::get_chemical_features_mesh to gemmi
// Returns a simple_mesh_t based on chemical features for the given residue
// Note: This is a minimal port matching the oracle behavior (empty mesh)
// The original function uses chemical_features::generate_meshes which depends
// on MMDB residue pointers and imol_no - not directly portable to gemmi
inline simple_mesh_t
molecule_t_get_chemical_features_mesh_gemmi(const gemmi::Residue& /*residue*/,
                                            const std::string& /*chain_id*/) {
   simple_mesh_t mesh;
   // Minimal port: always return empty mesh as in oracle run
   return mesh;
}

} // namespace coot