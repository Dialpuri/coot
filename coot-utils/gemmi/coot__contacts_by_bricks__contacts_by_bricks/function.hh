#pragma once

#include <gemmi/model.hpp>
#include <set>
#include <vector>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__contacts_by_bricks__set_lower_left_and_range/gemmi/function.hh"

namespace coot {
namespace contacts_by_bricks {

struct contacts_by_bricks {
    const gemmi::Atom* const* atoms;
    int n_atoms;
    float brick_size;
    float dist_nbc_max;
    bool only_between_different_residues_flag;
    float lower_left[3];
    int range[3];
    std::vector<std::vector<int>> atoms_in_bricks;
    std::vector<bool> fixed_flags;
    std::vector<std::vector<unsigned int>> thread_index_sets;
};

// No port for get_max_number_of_threads — inline translation
static unsigned int get_max_number_of_threads() {
    return 1u;
}

// No port for split_indices — inline translation
static void split_indices(
    std::vector<std::vector<unsigned int>>* indices,
    unsigned int n, unsigned int n_thread_sets) {
    indices->resize(n_thread_sets);
    unsigned int set_size = n / n_thread_sets;
    unsigned int remainder = n % n_thread_sets;
    unsigned int start = 0;
    for (unsigned int i = 0; i < n_thread_sets; i++) {
        unsigned int count = set_size + (i < remainder ? 1 : 0);
        for (unsigned int j = start; j < start + count; j++) {
            (*indices)[i].push_back(j);
        }
        start += count;
    }
}

// gemmi port of the contacts_by_bricks constructor
contacts_by_bricks contacts_by_bricks_gemmi(
    const gemmi::Atom* atoms_in[],
    int n_atoms_in,
    const std::set<unsigned int>& fixed_atom_indices)
{
    contacts_by_bricks result;
    result.atoms = atoms_in;
    result.n_atoms = n_atoms_in;
    result.brick_size = 20.0f;
    result.dist_nbc_max = 8.0f;
    result.only_between_different_residues_flag = false;

    result.range[0] = 0;
    result.range[1] = 0;
    result.range[2] = 0;
    set_lower_left_and_range_gemmi(atoms_in, n_atoms_in, result.brick_size, result.lower_left, result.range);

    int n_bricks = result.range[0] * result.range[1] * result.range[2];
    result.atoms_in_bricks.resize(n_bricks);

    // fill_the_bricks (inline gemmi translation)
    for (int i = 0; i < n_atoms_in; i++) {
        const gemmi::Atom* atom = atoms_in[i];
        if (atom) {
            int idx_3d[3];
            idx_3d[0] = static_cast<int>((atom->pos.x - result.lower_left[0]) / result.brick_size);
            idx_3d[1] = static_cast<int>((atom->pos.y - result.lower_left[1]) / result.brick_size);
            idx_3d[2] = static_cast<int>((atom->pos.z - result.lower_left[2]) / result.brick_size);
            int idx_1d = idx_3d[0] + idx_3d[1] * result.range[0] + idx_3d[2] * result.range[0] * result.range[1];
            result.atoms_in_bricks[idx_1d].push_back(i);
        }
    }

    result.fixed_flags.resize(n_atoms_in, false);
    for (std::set<unsigned int>::const_iterator it = fixed_atom_indices.begin();
         it != fixed_atom_indices.end(); ++it) {
        result.fixed_flags[*it] = true;
    }

    unsigned int n_threads = get_max_number_of_threads();
    unsigned int n_thread_sets = n_threads - 1;
    if (n_thread_sets < 1)
        n_thread_sets = 1;
    split_indices(&result.thread_index_sets, static_cast<unsigned int>(n_bricks), n_thread_sets);

    return result;
}

} // namespace contacts_by_bricks
} // namespace coot