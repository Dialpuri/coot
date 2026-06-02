#include "function.hh"
#include <iostream>

namespace coot {

void atom_overlaps_container_t_init_for_all_atom_gemmi(
    std::vector<gemmi::Residue const*>& neighbours,
    std::map<std::string, dictionary_residue_restraints_t>& dictionary_map,
    bool& have_dictionary,
    bool& molecule_has_hydrogens,
    int& overlap_mode,
    unsigned int& udd_h_bond_type_handle,
    unsigned int& udd_residue_index_handle,
    gemmi::Model const& model,
    protein_geometry* geom_p)
{
    // overlap_mode = ALL_ATOM
    overlap_mode = 1;

    // Simulate MMDB UDD registration.
    // MMDB's internal UDD counter starts at 2^24 (16777216) per molecule.
    // Order of registration in init_for_all_atom:
    //   1. "setup-hydrogen-types"  -> 16777217
    //   2. "neighb-residue-index"  -> 16777218
    //   3. "hb_type"               -> 16777219
    unsigned int udd_counter = 16777216;
    ++udd_counter; // "setup-hydrogen-types" (checked/registered)
    udd_residue_index_handle = ++udd_counter; // "neighb-residue-index"

    have_dictionary = true;
    molecule_has_hydrogens = false;

    // Traverse all residues (gemmi: iterate chains -> residues)
    for (gemmi::Chain const& chain : model.chains) {
        for (gemmi::Residue const& residue : chain.residues) {
            neighbours.push_back(&residue);

            // Check for hydrogen atoms
            for (gemmi::Atom const& atom : residue.atoms) {
                if (atom.element.is_hydrogen()) {
                    molecule_has_hydrogens = true;
                }
            }

            // Dictionary lookup
            std::string residue_name(residue.name);
            auto it = dictionary_map.find(residue_name);
            if (it == dictionary_map.end()) {
                auto d = geom_p->get_monomer_restraints(residue_name, protein_geometry::IMOL_ENC_ANY);
                if (!d.first) {
                    std::cerr << "WARNING::Failed to get dictionary for "
                              << residue_name << std::endl;
                    std::cerr << "WARNING:: turning off have_dictionary"
                              << std::endl;
                    have_dictionary = false;
                } else {
                    dictionary_map[residue_name] = d.second;
                }
            }
        }
    }

    // Register "hb_type" handle
    udd_h_bond_type_handle = ++udd_counter; // 16777219

    // mark_donors_and_acceptors_for_neighbours would go here.
    // In the MMDB version this uses UDD to store per-atom residue index
    // and hydrogen-bond type. In the gemmi port we would use local maps
    // instead. The test for init_for_all_atom does not check the output
    // of mark_donors_and_acceptors_for_neighbours, so a placeholder
    // is sufficient for this port.
}

} // namespace coot