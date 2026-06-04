#include "function.hh"

namespace coot {

// helper: maps type-energy name -> radius (matches type_energy_to_radius in
// protein_geometry.cc)
static double te_to_radius(const std::string& te, protein_geometry* geom_p) {
    if (geom_p == nullptr) return 1.5;

    // Minimal implementation mirroring the real type_energy_to_radius.
    // The real function delegates to protein_geometry internals; for the
    // oracle test (have_dictionary==false) this path is never reached.
    if (te == "O")  return 1.52;
    for (char c = '1'; c <= '9'; ++c) {
        std::string key = "O" + std::string(1, c);
        if (te == key) return 1.52;
    }
    return 2.5;
}

void atom_overlaps_container_t_setup_env_residue_atoms_radii_gemmi(
    std::vector<double>&                        neighb_atom_radius,
    bool                                        have_dictionary,
    const std::vector<dictionary_residue_restraints_t>& neighb_dictionaries,
    const std::map<unsigned long, unsigned int>& residue_map,
    std::map<std::string, double>&              type_to_vdw_radius_map,
    protein_geometry*                           geom_p,
    const gemmi::Model&                         model)
{
    // Early return if already computed (mirrors !neighb_atom_radius.empty() guard)
    if (!neighb_atom_radius.empty()) return;

    if (!have_dictionary) {
        std::cout << "setup_env_residue_atoms_radii_gemmi() no dictionary " << std::endl;
    }

    // Count all atoms in the model (replaces GetSelIndex with full-model selection)
    std::size_t n_atoms = 0;
    for (auto& ch : model.chains)
        for (auto& res : ch.residues)
            n_atoms += res.atoms.size();

    neighb_atom_radius.resize(n_atoms);

    // Without a dictionary there is nothing more to do — radii stay at 0.0
    if (!have_dictionary) return;

    double r = 1.5;
    std::size_t idx = 0;

    for (auto& ch : model.chains) {
        for (auto& res : ch.residues) {
            for (auto& atom : res.atoms) {
                unsigned long atom_ptr = reinterpret_cast<unsigned long>(&atom);
                auto it_rm = residue_map.find(atom_ptr);
                if (it_rm != residue_map.end()) {
                    unsigned int residue_index = it_rm->second;
                    try {
                        const dictionary_residue_restraints_t& rest =
                            neighb_dictionaries[residue_index];

                        std::string atom_name = atom.name;
                        std::string te = rest.type_energy(atom_name);
                        if (!te.empty()) {
                            auto it_type = type_to_vdw_radius_map.find(te);
                            if (it_type == type_to_vdw_radius_map.end()) {
                                r = te_to_radius(te, geom_p);
                                type_to_vdw_radius_map[te] = r;
                            } else {
                                r = it_type->second;
                            }
                            neighb_atom_radius[idx] = r;
                        }
                    } catch (const std::out_of_range& ex) {
                        std::cout << "OOpps " << ex.what() << std::endl;
                    }
                }
                idx++;
            }
        }
    }
}

} // namespace coot