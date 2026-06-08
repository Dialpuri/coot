#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <map>

namespace coot {

// Lightweight dict_torsion_restraint_t for gemmi port
struct dict_torsion_restraint_t {
    std::string name;
    std::string atom1;
    std::string atom2;
    std::string atom3;
    std::string atom4;
    double angle_;
    double esd_;
    int period_;

    dict_torsion_restraint_t(std::string n,
                             std::string a1, std::string a2,
                             std::string a3, std::string a4,
                             double ang, double esd, int per)
        : name(std::move(n)), atom1(std::move(a1)), atom2(std::move(a2)),
          atom3(std::move(a3)), atom4(std::move(a4)),
          angle_(ang), esd_(esd), period_(per) {}

    int periodicity() const { return period_; }
    double esd() const { return esd_; }
    std::string atom_id_1_4c() const { return atom1; }
    std::string atom_id_2_4c() const { return atom2; }
    std::string atom_id_3_4c() const { return atom3; }
    std::string atom_id_4_4c() const { return atom4; }
    double angle() const { return angle_; }
};

// Lightweight restraints_container_t for gemmi port
struct restraints_container_t {
    std::vector<std::vector<int>> restraints_vec;

    std::vector<bool> make_fixed_flags(int, int, int, int) {
        return {false, false, false, false};
    }

    void add(int type, int idx_1, int idx_2, int idx_3, int idx_4,
             std::vector<bool> fixed_flags,
             double angle, double esd, double weight, int period) {
        restraints_vec.push_back({type, idx_1, idx_2, idx_3, idx_4});
    }

    bool add_torsion_internal_gemmi(
        const dict_torsion_restraint_t &torsion_restraint,
        const std::vector<gemmi::Atom*> &res_atoms,
        const std::map<gemmi::Atom*, int> &atom_index_map,
        const double &torsion_restraint_weight);
};

}  // namespace coot
