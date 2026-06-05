#pragma once

#include <gemmi/model.hpp>
#include <set>
#include <vector>

namespace coot {
namespace contacts_by_bricks {

// ── Struct (copied from existing gemmi port of contacts_by_bricks) ──

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

static unsigned int get_max_number_of_threads() { return 1u; }

static void split_indices(
    std::vector<std::vector<unsigned int>>* indices,
    unsigned int n, unsigned int n_thread_sets)
{
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

// ── set_lower_left_and_range_gemmi ──

inline void set_lower_left_and_range_gemmi(
    gemmi::Atom const* atoms[],
    int n_atoms,
    float brick_size,
    float lower_left[3],
    int range[3])
{
    for (int i = 0; i < 3; i++)
        lower_left[i] = 9999.9f;

    unsigned int n = static_cast<unsigned int>(n_atoms);
    for (unsigned int i = 0; i < n; i++) {
        const gemmi::Atom* atom = atoms[i];
        float pos[3];
        pos[0] = atom->pos.x; pos[1] = atom->pos.y; pos[2] = atom->pos.z;
        for (int j = 0; j < 3; j++)
            if (pos[j] < lower_left[j])
                lower_left[j] = pos[j];
    }

    for (int i = 0; i < 3; i++)
        lower_left[i] -= 6.0f;

    float inv_brick_size = 1.0f / brick_size;
    for (unsigned int i = 0; i < n; i++) {
        const gemmi::Atom* atom = atoms[i];
        float pos[3];
        pos[0] = atom->pos.x; pos[1] = atom->pos.y; pos[2] = atom->pos.z;
        for (int j = 0; j < 3; j++) {
            float f = (pos[j] - lower_left[j]) * inv_brick_size;
            int brick_idx = static_cast<int>(f);
            if (brick_idx > range[j])
                range[j] = brick_idx;
        }
    }

    for (int i = 0; i < 3; i++)
        range[i] += 1;
}

// ── contacts_by_bricks_gemmi (constructor free function) ──

inline contacts_by_bricks contacts_by_bricks_gemmi(
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
    set_lower_left_and_range_gemmi(
        const_cast<gemmi::Atom const**>(atoms_in), n_atoms_in,
        result.brick_size, result.lower_left, result.range);

    int n_bricks = result.range[0] * result.range[1] * result.range[2];
    result.atoms_in_bricks.resize(n_bricks);

    for (int i = 0; i < n_atoms_in; i++) {
        const gemmi::Atom* atom = atoms_in[i];
        if (atom) {
            int idx_3d[3];
            idx_3d[0] = static_cast<int>((atom->pos.x - result.lower_left[0]) / result.brick_size);
            idx_3d[1] = static_cast<int>((atom->pos.y - result.lower_left[1]) / result.brick_size);
            idx_3d[2] = static_cast<int>((atom->pos.z - result.lower_left[2]) / result.brick_size);
            int idx_1d = idx_3d[0] + idx_3d[1] * result.range[0]
                       + idx_3d[2] * result.range[0] * result.range[1];
            result.atoms_in_bricks[idx_1d].push_back(i);
        }
    }

    result.fixed_flags.resize(n_atoms_in, false);
    for (auto it = fixed_atom_indices.begin(); it != fixed_atom_indices.end(); ++it) {
        result.fixed_flags[*it] = true;
    }

    unsigned int n_threads = get_max_number_of_threads();
    unsigned int n_thread_sets = n_threads - 1;
    if (n_thread_sets < 1)
        n_thread_sets = 1;
    split_indices(&result.thread_index_sets,
                  static_cast<unsigned int>(n_bricks), n_thread_sets);

    return result;
}

// ── find_the_contacts_between_bricks_simple_gemmi ──
// Scans neighbouring bricks and records atom-pair contacts within dist_nbc_max.
// atom_residue_ids[i] is a unique integer ID for the residue that atom i belongs to,
// used to implement the only_between_different_residues_flag check.

inline void find_the_contacts_between_bricks_simple_gemmi(
    const contacts_by_bricks& cbb,
    std::vector<std::set<unsigned int>>* vec,
    bool only_between_different_residues_flag,
    const std::vector<unsigned int>& atom_residue_ids)
{
    int n_bricks = cbb.range[0] * cbb.range[1] * cbb.range[2];
    int brick_index_max = n_bricks;
    float dist_max_sqrd = cbb.dist_nbc_max * cbb.dist_nbc_max;

    for (int ib = 0; ib < n_bricks; ib++) {
        const std::vector<int>& brick_base = cbb.atoms_in_bricks[ib];
        if (brick_base.size() > 0) {
            for (int iz = -1; iz < 2; iz++) {
                for (int iy = -1; iy < 2; iy++) {
                    for (int ix = -1; ix < 2; ix++) {
                        int ib_neighb = ib + ix + iy * cbb.range[0]
                                      + iz * cbb.range[0] * cbb.range[1];
                        if ((ib_neighb >= 0) && (ib_neighb != ib)) {
                            if (ib_neighb < brick_index_max) {
                                const std::vector<int>& brick_neighb =
                                    cbb.atoms_in_bricks[ib_neighb];
                                for (int base_idx : brick_base) {
                                    if (!cbb.fixed_flags[base_idx]) {
                                        const gemmi::Atom* at_1 = cbb.atoms[base_idx];
                                        for (int neighb_idx : brick_neighb) {
                                            if (only_between_different_residues_flag) {
                                                if (atom_residue_ids[static_cast<unsigned int>(neighb_idx)] ==
                                                    atom_residue_ids[static_cast<unsigned int>(base_idx)])
                                                    continue;
                                            }
                                            const gemmi::Atom* at_2 =
                                                cbb.atoms[neighb_idx];
                                            float d_x = at_1->pos.x - at_2->pos.x;
                                            float d_y = at_1->pos.y - at_2->pos.y;
                                            float d_z = at_1->pos.z - at_2->pos.z;
                                            float dd = d_x * d_x + d_y * d_y + d_z * d_z;
                                            if (dd < dist_max_sqrd) {
                                                vec->at(static_cast<unsigned int>(base_idx))
                                                    .insert(static_cast<unsigned int>(neighb_idx));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

} // namespace contacts_by_bricks
} // namespace coot