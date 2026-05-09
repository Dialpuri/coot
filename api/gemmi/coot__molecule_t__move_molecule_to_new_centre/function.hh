#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <utility>  // for std::pair

namespace coot {

// Helper struct to mimic clipper::Coord_orth behavior
struct Cartesian {
    float x_val, y_val, z_val;
    
    Cartesian() : x_val(0), y_val(0), z_val(0) {}
    Cartesian(float x, float y, float z) : x_val(x), y_val(y), z_val(z) {}
    
    float x() const { return x_val; }
    float y() const { return y_val; }
    float z() const { return z_val; }
    
    Cartesian operator-(const Cartesian& other) const {
        return Cartesian(x_val - other.x_val, y_val - other.y_val, z_val - other.z_val);
    }
    
    Cartesian operator+(const Cartesian& other) const {
        return Cartesian(x_val + other.x_val, y_val + other.y_val, z_val + other.z_val);
    }
};

// Helper: center of molecule (returns (is_valid, center))
// This mimics coot::centre_of_molecule behavior for gemmi::Structure
inline std::pair<bool, Cartesian> centre_of_molecule_gemmi(const gemmi::Structure& st) {
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    int n_atoms = 0;
    
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            for (const auto& residue : chain.residues) {
                for (const auto& atom : residue.atoms) {
                    sum_x += atom.pos.x;
                    sum_y += atom.pos.y;
                    sum_z += atom.pos.z;
                    n_atoms++;
                }
            }
        }
    }
    
    if (n_atoms == 0) {
        return std::make_pair(false, Cartesian(0, 0, 0));
    }
    
    return std::make_pair(true, Cartesian(sum_x / n_atoms, sum_y / n_atoms, sum_z / n_atoms));
}

// The ported function
inline int move_molecule_to_new_centre_gemmi(gemmi::Structure& st, const Cartesian& new_centre) {
    int status = 0;
    
    std::pair<bool, Cartesian> cm = centre_of_molecule_gemmi(st);
    if (cm.first) {
        Cartesian delta = new_centre - Cartesian(cm.second.x(), cm.second.y(), cm.second.z());
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& residue : chain.residues) {
                    for (auto& atom : residue.atoms) {
                        atom.pos.x += delta.x();
                        atom.pos.y += delta.y();
                        atom.pos.z += delta.z();
                    }
                }
            }
        }
    }
    
    return status;
}

} // namespace coot