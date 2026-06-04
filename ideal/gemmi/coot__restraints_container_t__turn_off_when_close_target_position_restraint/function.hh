#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <vector>
#include <atomic>
#include <thread>
#include <chrono>

// Restraint type constants (from coot::simple_restraint.hh)
constexpr int TARGET_POS_RESTRAINT = 6;

// ----------------------------------------------------------------
// simple_restraint_gemmi — minimal gemmi version of simple_restraint
// ----------------------------------------------------------------
struct simple_restraint_gemmi {
  int restraint_type;
  bool is_closed;
  int atom_index_1;
  gemmi::Vec3 atom_pull_target_pos;

  simple_restraint_gemmi()
    : restraint_type(0), is_closed(false), atom_index_1(-1),
      atom_pull_target_pos(0, 0, 0) {}

  void close() {
    is_closed = true;
  }
};

// ----------------------------------------------------------------
// restraints_container_t_gemmi — gemmi version with the one method we port
// ----------------------------------------------------------------
struct restraints_container_t_gemmi {
  std::vector<simple_restraint_gemmi> restraints_vec;
  std::vector<const gemmi::Atom*> atom;
  std::atomic<bool> restraints_lock;
  bool needs_reset;

  restraints_container_t_gemmi()
    : restraints_lock(false), needs_reset(false) {}

  bool turn_off_when_close_target_position_restraint_gemmi();
};

// ----------------------------------------------------------------
// turn_off_when_close_target_position_restraint_gemmi implementation
// ----------------------------------------------------------------
inline bool restraints_container_t_gemmi::turn_off_when_close_target_position_restraint_gemmi() {
  bool actioned = false;

  double close_dist = 0.6;

  bool unlocked = false;
  while (!restraints_lock.compare_exchange_weak(unlocked, true)) {
    std::this_thread::sleep_for(std::chrono::nanoseconds(100));
    unlocked = false;
  }

  for (auto &it : restraints_vec) {
    if (it.restraint_type == TARGET_POS_RESTRAINT) {
      if (it.is_closed) {
        if (static_cast<size_t>(it.atom_index_1) < atom.size()) {
          const gemmi::Atom *at = atom[static_cast<size_t>(it.atom_index_1)];
          if (at) {
            gemmi::Vec3 pos(at->pos.x, at->pos.y, at->pos.z);
            gemmi::Vec3 diff = pos - it.atom_pull_target_pos;
            double d = std::sqrt(diff.dot(diff));
            if (d < close_dist) {
              it.close();
              actioned = true;
            }
          }
        }
      }
    }
  }
  needs_reset = true;
  restraints_lock = false;
  return actioned;
}