#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_alt_confs/gemmi/function.hh"

namespace coot {
namespace primitive_chi_angles {

struct atom_name_quad {
    std::string names[4];
    atom_name_quad(const std::string& n1, const std::string& n2,
                   const std::string& n3, const std::string& n4) {
        names[0] = n1; names[1] = n2; names[2] = n3; names[3] = n4;
    }
    const std::string& atom_name(int i) const { return names[i]; }
};

struct atom_index_quad {
    int index1, index2, index3, index4;
    atom_index_quad(int i1, int i2, int i3, int i4)
        : index1(i1), index2(i2), index3(i3), index4(i4) {}
};

struct alt_confed_atom_index_quad {
    std::string alt_conf;
    std::vector<atom_index_quad> quad;
    alt_confed_atom_index_quad(const std::string& ac,
                               const std::vector<atom_index_quad>& q)
        : alt_conf(ac), quad(q) {}
};

inline std::vector<alt_confed_atom_index_quad>
get_quads_using_altconfs(const std::vector<atom_name_quad>& atom_name_quads,
                         const gemmi::Residue& residue) {
    std::vector<alt_confed_atom_index_quad> alt_v;
    
    const std::vector<std::string>& residue_alt_confs = 
        coot::util::get_residue_alt_confs_gemmi(residue);
    
    // Copy residue_alt_confs, removing empty strings
    std::vector<std::string> alt_confs;
    for (const auto& conf : residue_alt_confs) {
        if (!conf.empty()) {
            alt_confs.push_back(conf);
        }
    }
    
    // Handle the case where there are no alt confs (only empty string)
    if (alt_confs.empty()) {
        alt_confs.push_back("");
    }
    
    for (const auto& alt_conf : alt_confs) {
        std::vector<atom_index_quad> v;
        for (const auto& quad : atom_name_quads) {
            for (size_t i1 = 0; i1 < residue.atoms.size(); ++i1) {
                const auto& atom1 = residue.atoms[i1];
                if (atom1.name == quad.atom_name(0)) {
                    for (size_t i2 = 0; i2 < residue.atoms.size(); ++i2) {
                        const auto& atom2 = residue.atoms[i2];
                        if (atom2.name == quad.atom_name(1)) {
                            for (size_t i3 = 0; i3 < residue.atoms.size(); ++i3) {
                                const auto& atom3 = residue.atoms[i3];
                                if (atom3.name == quad.atom_name(2)) {
                                    for (size_t i4 = 0; i4 < residue.atoms.size(); ++i4) {
                                        const auto& atom4 = residue.atoms[i4];
                                        if (atom4.name == quad.atom_name(3)) {
                                            char ac1 = atom1.altloc;
                                            char ac2 = atom2.altloc;
                                            char ac3 = atom3.altloc;
                                            char ac4 = atom4.altloc;
                                            char target_ac = alt_conf.empty() ? ' ' : alt_conf[0];
                                            if ((ac4 == target_ac || ac4 == ' ') &&
                                                (ac3 == target_ac || ac3 == ' ') &&
                                                (ac2 == target_ac || ac2 == ' ') &&
                                                (ac1 == target_ac || ac1 == ' ')) {
                                                v.push_back(atom_index_quad(i1, i2, i3, i4));
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
        if (!v.empty()) {
            alt_v.push_back(alt_confed_atom_index_quad(alt_conf, v));
        }
    }
    
    return alt_v;
}

} // namespace primitive_chi_angles
} // namespace coot