#pragma once
#include <utility>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {

struct Cartesian {
    float x_val, y_val, z_val;
    
    Cartesian() : x_val(0), y_val(0), z_val(0) {}
    Cartesian(float x, float y, float z) : x_val(x), y_val(y), z_val(z) {}
    
    float x() const { return x_val; }
    float y() const { return y_val; }
    float z() const { return z_val; }
};

inline std::pair<bool, Cartesian> centre_of_molecule_gemmi(const gemmi::Structure& st) {
    float sum_x = 0.0f, sum_y = 0.0f, sum_z = 0.0f;
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
    
    float sf = 1.0f / static_cast<float>(n_atoms);
    return std::make_pair(true, Cartesian(sum_x * sf, sum_y * sf, sum_z * sf));
}

inline Cartesian
molecule_t_get_molecule_centre_gemmi(const gemmi::Structure& st) {
    std::pair<bool, Cartesian> centre = centre_of_molecule_gemmi(st);
    if (centre.first) {
        return Cartesian(centre.second.x(), centre.second.y(), centre.second.z());
    }
    return Cartesian(0.f, 0.f, 0.f);
}

} // namespace coot