#include "function.hh"
#include "ideal/refinement-results-mini-stats.hh"
#include <set>

namespace {

int find_atom_index_for_target_position(
    const coot::restraints_container_t& rc,
    const coot::atom_spec_t& atom_spec,
    bool& is_fixed)
{
    // Look up the atom in atom_map
    auto it = rc.atom_map.find(atom_spec);
    if (it != rc.atom_map.end()) {
        int idx = it->second;
        is_fixed = rc.fixed_check(idx);
        return idx;
    }
    return -1;
}

} // namespace

namespace coot {

void restraints_container_t::add_user_defined_target_position_restraint(
    int type, int atom_index, const atom_spec_t& atom_spec,
    const clipper::Coord_orth& pos, double weight)
{
    std::vector<int> restraint;
    restraint.push_back(type);
    restraint.push_back(atom_index);
    restraints_vec.push_back(restraint);
}

void add_extra_target_position_restraints_gemmi(
    restraints_container_t& rc,
    const extra_restraints_t& extra_restraints)
{
    for (const auto& pr : extra_restraints.target_position_restraints) {
        bool is_fixed = false;
        int atom_index =
            find_atom_index_for_target_position(rc, pr.atom_spec, is_fixed);

        if (!is_fixed && atom_index >= 0) {
            rc.add_user_defined_target_position_restraint(
                TARGET_POS_RESTRAINT,
                atom_index,
                pr.atom_spec,
                pr.pos,
                pr.weight);
        }
    }
}

} // namespace coot