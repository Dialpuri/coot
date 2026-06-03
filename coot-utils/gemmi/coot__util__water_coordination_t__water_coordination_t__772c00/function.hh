#pragma once

#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <gemmi/elem.hpp>
#include <vector>
#include <string>
#include <map>

namespace coot {
namespace util {

struct contact_atoms_info_t {
  // Central atom (the water/metal atom at the center)
  std::string central_chain;
  int central_res_seq;
  std::string central_res_name;
  std::string central_atom_name;

  // Contactors (atoms near the central atom)
  struct contactor_t {
    std::string chain;
    int res_seq;
    std::string res_name;
    std::string atom_name;
    double dist;
    gemmi::Element element;
  };
  std::vector<contactor_t> contactors;

  int n_contacts(double dist_max) const {
    int count = 0;
    for (const auto& c : contactors) {
      if (c.dist <= dist_max) count++;
    }
    return count;
  }

  bool has_metal_contact() const {
    for (const auto& c : contactors) {
      if (c.element.is_metal()) return true;
    }
    return false;
  }
};

class water_coordination_t {
public:
  std::vector<contact_atoms_info_t> contacts_;
  double radius_;
  bool do_metals_only_flag_;
  gemmi::Model* mol_;

  static water_coordination_t water_coordination_t_gemmi(gemmi::Structure* st,
                                                         double radius,
                                                         bool do_metals_only_flag) {
    water_coordination_t wc;
    if (!st || st->models.empty()) return wc;
    wc.init_internal(&st->models[0], st->cell, radius, do_metals_only_flag);
    return wc;
  }

  std::vector<contact_atoms_info_t> get_contacts() const {
    // Return only central atoms that have at least one contactor
    std::vector<contact_atoms_info_t> result;
    for (const auto& entry : contacts_) {
      if (!entry.contactors.empty()) {
        result.push_back(entry);
      }
    }
    return result;
  }

  std::vector<contact_atoms_info_t> metals() const {
    // Return central atoms that have metal contacts among their contactors
    std::vector<contact_atoms_info_t> result;
    for (const auto& entry : contacts_) {
      if (entry.has_metal_contact()) {
        result.push_back(entry);
      }
    }
    return result;
  }

  std::vector<contact_atoms_info_t> get_highly_coordinated_waters(int min_contacts, double dist_max) const {
    std::vector<contact_atoms_info_t> result;
    for (const auto& entry : contacts_) {
      if (entry.n_contacts(dist_max) >= min_contacts) {
        result.push_back(entry);
      }
    }
    return result;
  }

private:
  void init_internal(gemmi::Model* mol, const gemmi::UnitCell& cell, double radius, bool do_metals_only_flag) {
    if (!mol) return;

    mol_            = mol;
    radius_         = radius;
    do_metals_only_flag_ = do_metals_only_flag;

    // Use non-periodic cell (no symmetry) to match original MMDB behavior (use_symm=0)
    gemmi::UnitCell non_periodic;
    gemmi::NeighborSearch ns(*mol, non_periodic, radius + 0.5);
    ns.populate(false);

    // Collect central atoms: water (HOH) or metal atoms depending on mode
    struct central_atom_info {
      gemmi::CRA cra;
      gemmi::Position pos;
    };

    std::vector<central_atom_info> central_atoms;

    for (auto& chain : mol->chains) {
      for (auto& residue : chain.residues) {
        bool is_hoh = (residue.name == "HOH");

        if (!do_metals_only_flag) {
          // Normal mode: water atoms are central
          if (!is_hoh) continue;
          for (auto& atom : residue.atoms) {
            central_atoms.push_back({{&chain, &residue, &atom}, atom.pos});
          }
        } else {
          // Metals-only mode: metal atoms are central
          for (auto& atom : residue.atoms) {
            if (atom.element.is_metal()) {
              central_atoms.push_back({{&chain, &residue, &atom}, atom.pos});
            }
          }
        }
      }
    }

    // For each central atom, create an entry and find contacts
    for (auto& central : central_atoms) {
      contact_atoms_info_t entry;
      entry.central_chain  = central.cra.chain->name;
      entry.central_res_seq = central.cra.residue->seqid.num.value;
      entry.central_res_name = central.cra.residue->name;
      entry.central_atom_name = central.cra.atom->name;

      auto marks = ns.find_atoms(central.pos, '\0', 0.5, radius);

      for (auto* m : marks) {
        gemmi::CRA cra = m->to_cra(*mol);
        if (!cra.chain || !cra.residue || !cra.atom) continue;

        // Skip water residues themselves
        if (cra.residue->name == "HOH") continue;

        // Skip hydrogen atoms
        if (cra.atom->element.is_hydrogen()) continue;

        // Compute distance
        double dist = (cra.atom->pos - central.pos).length();
        if (dist > radius) continue;

        contact_atoms_info_t::contactor_t c{
          cra.chain->name,
          cra.residue->seqid.num.value,
          cra.residue->name,
          cra.atom->name,
          dist,
          cra.atom->element
        };

        entry.contactors.push_back(c);
      }

      contacts_.push_back(entry);
    }
  }
};

} // namespace util
} // namespace coot