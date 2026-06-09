#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <vector>
#include <string>
#include <cmath>

namespace coot {

// Restraint type enum (from coot source)
enum restraint_type_t {
   ANGLE_RESTRAINT = 0,
   BOND_RESTRAINT,
   PLANAR_RESTRAINT,
   TRANS_PEPTIDE_RESTRAINT,
   NON_BOND_BUMP_RESTRAINT,
   NON_BOND_CONTACT_RESTRAINT,
   HYDROGEN_BOND_RESTRAINT,
   CHIRAL_RESTRAINT,
   NON_BOND_VDW_RESTRAINT,
   NON_BOND_1_4_RESTRAINT,
   DISULFIDE_BOND_RESTRAINT,
   RESTRAINT_TYPE_END
};

// Simple restraint struct with is_closed field
struct simple_restraint_t {
   restraint_type_t rest_type;
   int atom_index_1;
   int atom_index_2;
   int atom_index_3;
   int atom_index_4;
   bool is_closed;
   float target;
   float sigma;
   float torsion_restraint_weight;
   int periodicity;
   std::vector<bool> fixed_atom_flags;

   simple_restraint_t(restraint_type_t rt, int a1, int a2, int a3, int a4,
                      bool closed, float t, float s, float w, int p)
      : rest_type(rt), atom_index_1(a1), atom_index_2(a2), atom_index_3(a3), atom_index_4(a4),
        is_closed(closed), target(t), sigma(s), torsion_restraint_weight(w), periodicity(p) {}
};

class restraints_container_t {
public:
   std::vector<simple_restraint_t> restraints_vec;
   std::vector<gemmi::CRA> atom_restraints_vec;

   unsigned int size() const { return restraints_vec.size(); }

   void post_add_new_restraint() {
      // Hook for post-processing; empty in gemmi port
   }

   void add_trans_peptide_restraint_gemmi(
       gemmi::Chain* first_chain, gemmi::Residue* first,
       gemmi::Chain* second_chain, gemmi::Residue* second);

private:
   // Compute torsion angle (dihedral) from 4 positions
   static double torsion_angle(const gemmi::Vec3& a, const gemmi::Vec3& b,
                                const gemmi::Vec3& c, const gemmi::Vec3& d) {
      gemmi::Vec3 b1 = b - a;
      gemmi::Vec3 b2 = c - b;
      gemmi::Vec3 b3 = d - c;

      // Manual cross product
      auto cross = [](const gemmi::Vec3& u, const gemmi::Vec3& v) -> gemmi::Vec3 {
          return gemmi::Vec3(
              u.y * v.z - u.z * v.y,
              u.z * v.x - u.x * v.z,
              u.x * v.y - u.y * v.x
          );
      };

      gemmi::Vec3 n1 = cross(b1, b2);
      gemmi::Vec3 n2 = cross(b2, b3);

      double n1_len = n1.length();
      double n2_len = n2.length();

      if (n1_len < 1e-10 || n2_len < 1e-10) return 0.0;

      n1 = gemmi::Vec3(n1.x / n1_len, n1.y / n1_len, n1.z / n1_len);
      n2 = gemmi::Vec3(n2.x / n2_len, n2.y / n2_len, n2.z / n2_len);

      double b2_len = b2.length();
      if (b2_len < 1e-10) return 0.0;
      b2 = gemmi::Vec3(b2.x / b2_len, b2.y / b2_len, b2.z / b2_len);

      gemmi::Vec3 m1 = cross(n1, b2);

      double x = n1.dot(n2);
      double y = m1.dot(n2);

      return std::atan2(y, x);
   }
};

inline void restraints_container_t::add_trans_peptide_restraint_gemmi(
    gemmi::Chain* first_chain, gemmi::Residue* first,
    gemmi::Chain* second_chain, gemmi::Residue* second) {

  // Does the restraint already exist?
  bool new_restraint_needed = true;

  for (auto& r : restraints_vec) {
    if (r.rest_type == TRANS_PEPTIDE_RESTRAINT) {
      gemmi::CRA a1 = atom_restraints_vec[r.atom_index_1];
      gemmi::CRA a2 = atom_restraints_vec[r.atom_index_2];
      gemmi::CRA a3 = atom_restraints_vec[r.atom_index_3];
      gemmi::CRA a4 = atom_restraints_vec[r.atom_index_4];

      gemmi::Residue* r_11 = a1.residue;
      gemmi::Residue* r_12 = a2.residue;
      gemmi::Residue* r_21 = a3.residue;
      gemmi::Residue* r_22 = a4.residue;

      if (r_11 == first) {
        if (r_12 == first) {
          if (r_21 == second) {
            if (r_22 == second) {
              r.is_closed = false;
              new_restraint_needed = false;
            }
          }
        }
      }
    }
  }

  if (new_restraint_needed) {
    // Find CA, C in first residue and N, CA in second residue
    const gemmi::Atom* ca1 = nullptr;
    const gemmi::Atom* c1  = nullptr;
    const gemmi::Atom* n2  = nullptr;
    const gemmi::Atom* ca2 = nullptr;

    for (const auto& atom : first->atoms) {
      if (atom.name == "CA") ca1 = &atom;
      if (atom.name == "C")  c1  = &atom;
    }
    for (const auto& atom : second->atoms) {
      if (atom.name == "N")  n2  = &atom;
      if (atom.name == "CA") ca2 = &atom;
    }

    if (ca1 && c1 && n2 && ca2) {
      double cn_dist = (c1->pos - n2->pos).length();
      bool add_restraint = (cn_dist < 2.0);

      if (add_restraint) {
        double omega = torsion_angle(ca1->pos, c1->pos, n2->pos, ca2->pos);
        bool approximately_trans = (std::abs(omega) > M_PI / 2.0);

        if (approximately_trans) {
          int idx1 = static_cast<int>(atom_restraints_vec.size());
          atom_restraints_vec.push_back(gemmi::CRA{first_chain, first, const_cast<gemmi::Atom*>(ca1)});

          int idx2 = static_cast<int>(atom_restraints_vec.size());
          atom_restraints_vec.push_back(gemmi::CRA{first_chain, first, const_cast<gemmi::Atom*>(c1)});

          int idx3 = static_cast<int>(atom_restraints_vec.size());
          atom_restraints_vec.push_back(gemmi::CRA{second_chain, second, const_cast<gemmi::Atom*>(n2)});

          int idx4 = static_cast<int>(atom_restraints_vec.size());
          atom_restraints_vec.push_back(gemmi::CRA{second_chain, second, const_cast<gemmi::Atom*>(ca2)});

          std::vector<bool> fixed_flags = {false, false, false, false};
          restraints_vec.push_back(
              simple_restraint_t(TRANS_PEPTIDE_RESTRAINT,
                                 idx1, idx2, idx3, idx4,
                                 true,      // is_closed
                                 static_cast<float>(omega),
                                 100.0f,    // sigma
                                 1.0f,      // torsion_restraint_weight
                                 1));       // periodicity
        }
      }
    }
    post_add_new_restraint();
  }
}

} // namespace coot
