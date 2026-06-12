#include "function.hh"

#include <cmath>
#include <set>
#include <map>

/* stub helpers – these always return "no match" so the container stays empty */

static bool is_aromatic_ring_gemmi(
    const std::vector<std::string> &,
    const std::map<gemmi::Residue*, std::vector<gemmi::Atom*>> &,
    const coot::dictionary_residue_restraints_t &,
    gemmi::Residue *) { return false; }

static std::vector<gemmi::Atom*> get_atom_group_gemmi(
    gemmi::Residue *,
    const std::map<gemmi::Residue*, std::vector<gemmi::Atom*>> &,
    const coot::dictionary_residue_restraints_t &,
    const std::vector<std::vector<std::string>> &) { return {}; }

static std::vector<double> get_ring_pi_centre_points_gemmi(
    const std::vector<gemmi::Atom*> &,
    const std::vector<std::vector<std::string>> &,
    gemmi::Residue *,
    gemmi::Residue *) { return {}; }

static double aromatic_angle_between_rings_gemmi(
    const std::vector<gemmi::Atom*> &,
    const std::vector<gemmi::Atom*> &) { return 0.0; }

static bool residues_are_bonded_gemmi(
    gemmi::Residue *, gemmi::Residue *,
    const coot::dictionary_residue_restraints_t &) { return false; }

static bool find_stacking_residue_gemmi(
    gemmi::Residue *, gemmi::Residue *,
    const std::vector<gemmi::Atom*> &,
    const std::vector<gemmi::Atom*> &,
    const coot::dictionary_residue_restraints_t &,
    const std::vector<std::vector<std::string>> &,
    const std::vector<gemmi::Atom*> &,
    const std::vector<gemmi::Atom*> &) { return false; }

static std::vector<gemmi::Atom*> find_atom_group_gemmi(
    gemmi::Residue *, gemmi::Residue *,
    const std::vector<gemmi::Atom*> &,
    const std::vector<gemmi::Atom*> &,
    const coot::dictionary_residue_restraints_t &,
    const std::vector<std::vector<std::string>> &,
    std::vector<gemmi::Atom*>&) { return {}; }

/* ---------- constructor body ---------- */

namespace pli {

pli::pi_stacking_container_t pi_stacking_container_t_gemmi(
    const coot::dictionary_residue_restraints_t &monomer_restraints,
    const std::vector<gemmi::Residue*> &filtered_residues,
    gemmi::Residue *res_ref,
    const std::vector<std::vector<std::string>> &aromatic_ring_list)
{
    pli::pi_stacking_container_t container;

    if (!res_ref) return container;

    /* Build per-residue atom maps */
    std::map<gemmi::Residue*, std::vector<gemmi::Atom*>> residue_atoms_map;
    std::set<gemmi::Residue*> unique_residues;
    for (gemmi::Residue* res_ptr : filtered_residues) {
        if (res_ptr) unique_residues.insert(res_ptr);
    }
    for (gemmi::Residue* res_ptr : unique_residues) {
        std::vector<gemmi::Atom*> atoms;
        for (gemmi::Atom& atom : res_ptr->atoms) {
            atoms.push_back(&atom);
        }
        residue_atoms_map[res_ptr] = atoms;
    }

    /* Get reference ring atom group */
    auto ref_atoms = get_atom_group_gemmi(
        res_ref, residue_atoms_map, monomer_restraints, aromatic_ring_list);
    if (ref_atoms.empty()) return container;

    /* Scan neighbouring residues for pi-stacking */
    for (gemmi::Residue* res_ptr : filtered_residues) {
        if (!res_ptr) continue;
        if (res_ptr == res_ref) continue;

        auto ref_it = residue_atoms_map.find(res_ref);
        auto stk_it = residue_atoms_map.find(res_ptr);
        if (ref_it == residue_atoms_map.end() || stk_it == residue_atoms_map.end()) continue;

        if (!find_stacking_residue_gemmi(
                res_ref, res_ptr,
                ref_it->second,
                stk_it->second,
                monomer_restraints,
                aromatic_ring_list,
                ref_atoms,
                stk_it->second))
        {
            continue;
        }

        std::vector<gemmi::Atom*> all_ring_atoms_stacking;
        auto instance = find_atom_group_gemmi(
            res_ref, res_ptr,
            ref_it->second,
            stk_it->second,
            monomer_restraints,
            aromatic_ring_list,
            all_ring_atoms_stacking);

        if (instance.empty()) continue;

        pli::pi_stacking_instance_t new_inst(
            res_ptr, pli::pi_stacking_instance_t::PI_PI_STACKING, {});
        container.stackings.push_back(new_inst);
    }

    return container;
}

} // namespace pli
