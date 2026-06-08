#pragma once

#include <vector>
#include <cstddef>
#include <cstdint>
#include <cmath>

namespace coot {

struct restraints_container_t_gemmi {
  unsigned int n_threads = 0;
  std::vector<std::vector<std::size_t>> restraints_indices;
  std::vector<std::vector<double>> df_by_thread_results;
  std::vector<std::vector<std::size_t>> df_by_thread_atom_indices;

  // Configurable for testing - matches what real restraints_container_t returns
  unsigned int restraint_count_ = 0;
  unsigned int n_vars_ = 0;
  unsigned int atom_count_ = 0;

  unsigned int size() const { return restraint_count_; }
  unsigned int n_variables() const { return n_vars_; }
  unsigned int get_n_atoms() const { return atom_count_; }

  void make_df_restraints_indices_gemmi();
};

inline void restraints_container_t_gemmi::make_df_restraints_indices_gemmi() {
   unsigned int n_r_s = n_threads;
   unsigned int restraints_size = size();

   if (n_r_s == 0) n_r_s = 1;

   restraints_indices.clear();
   restraints_indices.resize(n_r_s);
   if (df_by_thread_results.size() > 0) {
      df_by_thread_results.clear();
   }
   df_by_thread_results.resize(n_r_s);

   int r_reserve_size = std::lround(static_cast<float>(restraints_size)/static_cast<float>(n_r_s)) + 2;

   for (std::size_t n=0; n<restraints_indices.size(); n++) {
      restraints_indices[n].reserve(r_reserve_size);
   }

   unsigned int i_thread = 0;
   for (unsigned int ir=0; ir<restraints_size; ir++) {
      restraints_indices[i_thread].push_back(ir);
      ++i_thread;
      if (i_thread==n_r_s) i_thread=0;
   }

   unsigned int n_var = n_variables();
   for (std::size_t ii=0; ii<n_r_s; ii++)
      df_by_thread_results[ii] = std::vector<double>(n_var, 0);

   df_by_thread_atom_indices.clear();
   df_by_thread_atom_indices.resize(n_r_s);
   i_thread = 0;
   unsigned int n = get_n_atoms();
   for (unsigned int ir=0; ir<n; ir++) {
      df_by_thread_atom_indices[i_thread].push_back(ir);
      ++i_thread;
      if (i_thread==n_r_s) i_thread=0;
   }
}

} // namespace coot
