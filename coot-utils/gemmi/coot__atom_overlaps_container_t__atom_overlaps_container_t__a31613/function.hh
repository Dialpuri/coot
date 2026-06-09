#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <map>
#include <string>

#include "geometry/protein-geometry.hh"

namespace coot {

// gemmi-compatible overlap record (replaces MMDB atom_overlap_t with mmdb::Atom*)
struct atom_overlap_t_gemmi {
    double distance = 0.0;
    double weight = 0.0;
    gemmi::Atom const* atom1 = nullptr;
    gemmi::Atom const* atom2 = nullptr;
    int overlap_type = 0;
};

// Result struct mirroring atom_overlaps_container_t members (gemmi version)
struct atom_overlaps_container_result_t {
    gemmi::Residue const* res_central = nullptr;
    std::vector<gemmi::Residue const*> neighbours;
    int imol_enc = 0;
    double clash_spike_length = 0.0;
    double probe_radius = 0.0;
    const protein_geometry* geom_p = nullptr;

    std::vector<atom_overlap_t_gemmi> overlaps;
    int overlap_mode = 0;
    bool have_dictionary = false;
    bool molecule_has_hydrogens = false;
    std::vector<dictionary_residue_restraints_t> neighb_dictionaries;
    dictionary_residue_restraints_t central_residue_dictionary;
    int udd_h_bond_type_handle = 0;
    int udd_residue_index_handle = -1;
};

// gemmi port of the atom_overlaps_container_t constructor
// Signature mirrors the original but uses gemmi types
//
// Original:
//   coot::atom_overlaps_container_t::atom_overlaps_container_t(
//       mmdb::Residue *res_central_in,
//       const std::vector<mmdb::Residue *> &neighbours_in,
//       mmdb::Manager *mol_in,
//       int imol_enc_in,
//       const protein_geometry *geom_p_in,
//       double clash_spike_length_in,
//       double probe_radius_in)
//
// gemmi:
inline atom_overlaps_container_result_t atom_overlaps_container_t_gemmi(
    gemmi::Residue const* res_central_in,
    std::vector<gemmi::Residue const*> const& neighbours_in,
    gemmi::Model const& model,
    int imol_enc_in,
    const protein_geometry* geom_p_in,
    double clash_spike_length_in,
    double probe_radius_in)
{
    atom_overlaps_container_result_t result;

    // Constructor member initialisation (mirrors original)
    result.probe_radius        = probe_radius_in;
    result.geom_p              = geom_p_in;
    result.res_central         = res_central_in;
    result.neighbours          = neighbours_in;
    result.imol_enc            = imol_enc_in;
    result.clash_spike_length  = clash_spike_length_in;

    // --- init() body (translated to gemmi) ---

    // CENTRAL_RESIDUE constant from atom-overlaps.hh
    static const int CENTRAL_RESIDUE = 0;
    result.overlap_mode           = CENTRAL_RESIDUE;
    result.udd_residue_index_handle = -1; // unset
    result.have_dictionary        = false;
    result.molecule_has_hydrogens = false;

    if (res_central_in && geom_p_in) {
        // gemmi: residue.name is ResidueId, .name field is the 3-letter code
        std::string cres_name = res_central_in->name;

        auto d = geom_p_in->get_monomer_restraints(cres_name, imol_enc_in);

        if (!d.first) {
            // Warning path — original prints to std::cout
        } else {
            result.central_residue_dictionary = d.second;
            result.neighb_dictionaries.resize(neighbours_in.size());
            result.have_dictionary = true;

            for (unsigned int i = 0; i < neighbours_in.size(); i++) {
                std::string residue_name = neighbours_in[i]->name;
                d = geom_p_in->get_monomer_restraints(residue_name, imol_enc_in);
                if (!d.first) {
                    result.have_dictionary = false;
                    break;
                }
                result.neighb_dictionaries[i] = d.second;
            }
        }
        // NOTE: fill_ligand_atom_neighbour_map() and mark_donors_and_acceptors()
        // are called by the original when have_dictionary is true, but they
        // populate internal maps that don't affect the `overlaps` vector.
        // The overlaps vector is populated by check_overlaps(), not init().
    }

    return result;
}

} // namespace coot
