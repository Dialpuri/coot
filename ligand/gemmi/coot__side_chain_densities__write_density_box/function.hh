#pragma once

#include <gemmi/model.hpp>
#include <string>

namespace coot::side_chain_densities {

/// Gemmi port of coot::side_chain_densities::write_density_box.
/// Writes a 3D density box to a .tab file under side-chain-data/<res>/<rotamer>/<id>-<chain>-<seq>.tab
///
/// @param density_box  raw pointer to density values (caller-owned)
/// @param n_steps      half-width of grid in each dimension (box size = 2*n_steps+1)
/// @param cra          gemmi::CRA with at least cra.residue non-null; cra.chain used for chain name
/// @param id           identifier string (e.g. "test")
/// @param rotamer_name rotamer name (e.g. "m120"); empty string skips file write
void write_density_box_gemmi(
    float* density_box,
    int    n_steps,
    gemmi::CRA cra,
    const std::string& id,
    const std::string& rotamer_name);

} // namespace coot::side_chain_densities