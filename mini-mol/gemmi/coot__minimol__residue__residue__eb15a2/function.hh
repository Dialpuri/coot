#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {
namespace minimol {

class atom {
public:
    std::string altLoc;
    float occupancy;
    float temperature_factor;
    float x, y, z;
    std::string name;
    std::string element;
    int int_user_data;

    atom()
        : x(0), y(0), z(0), occupancy(0), temperature_factor(0),
          int_user_data(-1) {}
    atom(std::string atom_name, std::string ele, float x_c, float y_c, float z_c,
         const std::string &altloc, float occupancy_v, float dbf)
        : name(atom_name), element(ele), x(x_c), y(y_c), z(z_c),
          altLoc(altloc), occupancy(occupancy_v), temperature_factor(dbf),
          int_user_data(-1) {}
};

class residue {
public:
    int seqnum;
    std::string ins_code;
    std::string name;
    std::vector<atom> atoms;

    residue() : seqnum(-2147483647), ins_code(""), name("") {}

    explicit residue(int i) : seqnum(i), ins_code(""), name("") {}

    residue(int i, const std::string &resname)
        : seqnum(i), ins_code(""), name(resname) {}

    residue(const gemmi::Residue *residue_p);

    residue(const gemmi::Residue *residue_p,
            const std::vector<std::string> &keep_only_these_atoms);

    void addatom(std::string atom_name, std::string element,
                 float x_c, float y_c, float z_c, const std::string &altloc,
                 float bf, float occupancy_v) {
        atom a(atom_name, element, x_c, y_c, z_c, altloc, occupancy_v, bf);
        atoms.push_back(a);
    }

    unsigned int n_atoms() const { return atoms.size(); }

    const atom &operator[](int i) const { return atoms[i]; }
    atom &operator[](int i) { return atoms[i]; }
};

inline std::string normalize_ins_code(char c) {
    return (c == ' ') ? std::string("") : std::string(1, c);
}

inline std::string pad_atom_name_pdb(const std::string &name) {
    if (name.size() >= 4) return name;
    return std::string(1, ' ') + name +
           std::string(4 - 1 - name.size(), ' ');
}

// gemmi-based constructor from gemmi::Residue*
inline residue::residue(const gemmi::Residue *residue_p)
    : seqnum(residue_p->seqid.num.value),
      ins_code(normalize_ins_code(residue_p->seqid.icode)),
      name(residue_p->name) {
    for (const auto &g_atom : residue_p->atoms) {
        std::string atom_name = pad_atom_name_pdb(g_atom.name);
        addatom(atom_name, g_atom.element.name(), g_atom.pos.x, g_atom.pos.y,
                g_atom.pos.z, std::string(1, g_atom.altloc), g_atom.b_iso,
                g_atom.occ);
    }
}

// gemmi-based constructor from gemmi::Residue* with keep list
inline residue::residue(const gemmi::Residue *residue_p,
                        const std::vector<std::string> &keep_only_these_atoms)
    : seqnum(residue_p->seqid.num.value),
      ins_code(normalize_ins_code(residue_p->seqid.icode)),
      name(residue_p->name) {
    for (const auto &g_atom : residue_p->atoms) {
        std::string atom_name = pad_atom_name_pdb(g_atom.name);
        bool add_it = false;
        for (const auto &keep : keep_only_these_atoms) {
            if (atom_name == keep) {
                add_it = true;
                break;
            }
        }
        if (add_it) {
            addatom(atom_name, g_atom.element.name(), g_atom.pos.x,
                    g_atom.pos.y, g_atom.pos.z,
                    std::string(1, g_atom.altloc), g_atom.b_iso, g_atom.occ);
        }
    }
}

} // namespace minimol
} // namespace coot

// Factory functions named residue_gemmi to match the port naming convention
namespace coot {
namespace minimol {

inline residue residue_gemmi() { return residue(); }

inline residue residue_gemmi(int i) { return residue(i); }

inline residue residue_gemmi(int i, const std::string &resname) {
    return residue(i, resname);
}

inline residue residue_gemmi(const gemmi::Residue *residue_p) {
    return residue(residue_p);
}

inline residue residue_gemmi(const gemmi::Residue *residue_p,
                             const std::vector<std::string> &keep_only_these_atoms) {
    return residue(residue_p, keep_only_these_atoms);
}

} // namespace minimol
} // namespace coot
