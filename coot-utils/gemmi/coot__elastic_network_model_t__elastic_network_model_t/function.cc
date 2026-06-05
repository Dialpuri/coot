#include "function.hh"

#include <cmath>
#include <cstdlib>
#include <map>

namespace coot {

elastic_network_model_t_gemmi::elastic_network_model_t_gemmi(
    const std::vector<const gemmi::Atom*>& selected_atoms,
    double min_dist,
    double max_dist,
    unsigned int max_n_distances) {

  if (selected_atoms.empty() || max_n_distances == 0)
    return;

  // Find all contacts within [min_dist, max_dist] among selected atoms
  std::vector<std::pair<int, int>> contacts;
  int n = static_cast<int>(selected_atoms.size());
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      const gemmi::Atom* a1 = selected_atoms[i];
      const gemmi::Atom* a2 = selected_atoms[j];
      if (!a1 || !a2) continue;
      double dx = a1->pos.x - a2->pos.x;
      double dy = a1->pos.y - a2->pos.y;
      double dz = a1->pos.z - a2->pos.z;
      double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (dist >= min_dist && dist <= max_dist) {
        contacts.push_back(std::make_pair(i, j));
      }
    }
  }

  int n_contacts = static_cast<int>(contacts.size());
  if (n_contacts == 0) return;

  if (n_contacts <= static_cast<int>(max_n_distances)) {
    // Add all contacts
    for (int i = 0; i < n_contacts; i++) {
      const gemmi::Atom* at1 = selected_atoms[contacts[i].first];
      const gemmi::Atom* at2 = selected_atoms[contacts[i].second];
      elastic_network_item_t_gemmi item(at1, at2, 0.1f);
      d.push_back(item);
    }
  } else {
    // Randomly sample max_n_distances unique contacts
    double inv_rand_max = 1.0 / double(RAND_MAX);
    std::map<int, elastic_network_item_t_gemmi> contact_indices;
    while (static_cast<unsigned int>(contact_indices.size()) < max_n_distances) {
      int random_index = static_cast<int>(double(std::rand()) * inv_rand_max * double(n_contacts));
      random_index = std::abs(random_index) % n_contacts;
      if (contact_indices.find(random_index) == contact_indices.end()) {
        const gemmi::Atom* at1 = selected_atoms[contacts[random_index].first];
        const gemmi::Atom* at2 = selected_atoms[contacts[random_index].second];
        contact_indices[random_index] = elastic_network_item_t_gemmi(at1, at2, 0.1f);
      }
    }
    for (auto it = contact_indices.begin(); it != contact_indices.end(); ++it) {
      d.push_back(it->second);
    }
  }
}

} // namespace coot