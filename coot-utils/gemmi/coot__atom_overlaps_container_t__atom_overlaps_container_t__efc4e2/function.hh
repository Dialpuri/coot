#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <string>
#include <vector>

#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_overlaps_container_t__init_for_all_atom/gemmi/function.hh"

namespace coot {

/* gemmi-port of the atom_overlaps_container_t constructor.
 *
 * Mirrors the original behaviour:
 *   - sets geom_p, probe_radius, ignore_water_contacts_flag
 *   - hardcodes clash_spike_length to 0.5 (same as original)
 *   - calls init_for_all_atom (via its gemmi port)
 *
 * Returns a small result struct containing the gemmi-compatible members
 * that the original container would have populated.
 */
struct atom_overlaps_container_result_t {
    std::vector<gemmi::Residue const*> neighbours;
    std::map<std::string, dictionary_residue_restraints_t> dictionary_map;
    bool have_dictionary = false;
    bool molecule_has_hydrogens = false;
    int overlap_mode = 0;
    unsigned int udd_h_bond_type_handle = 0;
    unsigned int udd_residue_index_handle = 0;
    bool ignore_water_contacts_flag = false;
    double clash_spike_length = 0.0;
    double probe_radius = 0.0;
};

inline atom_overlaps_container_result_t atom_overlaps_container_t_gemmi(
    gemmi::Model const& model,
    const protein_geometry* geom_p_in,
    bool ignore_water_contacts_flag_in,
    double /*clash_spike_length_in*/,
    double probe_radius_in)
{
    atom_overlaps_container_result_t result;

    // Mirror the original constructor's member initialisation
    result.ignore_water_contacts_flag = ignore_water_contacts_flag_in;
    result.clash_spike_length = 0.5; // hardcoded in original, input ignored
    result.probe_radius = probe_radius_in;

    // Delegate to the already-portd init_for_all_atom
    atom_overlaps_container_t_init_for_all_atom_gemmi(
        result.neighbours,
        result.dictionary_map,
        result.have_dictionary,
        result.molecule_has_hydrogens,
        result.overlap_mode,
        result.udd_h_bond_type_handle,
        result.udd_residue_index_handle,
        model,
        const_cast<protein_geometry*>(geom_p_in));

    return result;
}

} // namespace coot