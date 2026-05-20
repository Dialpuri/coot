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
    
    atom() : x(0), y(0), z(0), occupancy(0), temperature_factor(0), int_user_data(-1) {}
    atom(std::string atom_name, std::string ele, float x, float y, float z, const std::string &altloc, float occupancy, float dbf)
        : name(atom_name), element(ele), x(x), y(y), z(z), altLoc(altloc), occupancy(occupancy), temperature_factor(dbf), int_user_data(-1) {}
};

class residue {
public:
    explicit residue(int i) : ins_code(""), name("") { seqnum = i; }
    residue(int i, const std::string& resname) : ins_code(""), name(resname) {
        seqnum = i;
    }
    residue() { seqnum = -2147483647; }
    residue(const gemmi::Residue* residue_p,
            const std::vector<std::string>& keep_only_these_atoms);
    residue(const gemmi::Residue& residue_ref,
            const std::vector<std::string>& keep_only_these_atoms);
    
    int seqnum;
    std::string ins_code;
    std::string name;
    std::vector<atom> atoms;
    
    void addatom(std::string atom_name, std::string element,
                 float x, float y, float z, const std::string& altloc, float bf, float occupancy);
    unsigned int n_atoms() const { return atoms.size(); }
    const atom& operator[](int i) const { return atoms[i]; }
    atom& operator[](int i) { return atoms[i]; }
};

// Helper to convert gemmi's ' ' (space) to "" for insertion code
inline std::string normalize_ins_code(char c) {
    return (c == ' ') ? std::string("") : std::string(1, c);
}

// Helper to left-pad atom name to 4 characters (PDB format)
// PDB atom name is 4 chars, left-aligned. gemmi stores "CA" without padding.
// We need " CA " (space + CA + space) but gemmi gives us "CA".
// For 2-char names, pad with 1 leading space + 1 trailing space: " CA "
inline std::string pad_atom_name_pdb(const std::string& name) {
    if (name.size() >= 4) return name;
    // PDB format: columns 13-16, atom name is left-aligned
    // "CA" in PDB is stored as " CA " (space at col 13, C at 14, A at 15, space at 16)
    return std::string(1, ' ') + name + std::string(4 - 1 - name.size(), ' ');
}

// Inline implementation
inline residue::residue(const gemmi::Residue* residue_p,
                        const std::vector<std::string>& keep_only_these_atoms) {
    seqnum = residue_p->seqid.num.value;
    ins_code = normalize_ins_code(residue_p->seqid.icode);
    name = residue_p->name;
    
    for (const auto& atom : residue_p->atoms) {
        // Pad atom name to 4 characters in PDB format (1 leading space + name + trailing spaces)
        std::string atom_name = pad_atom_name_pdb(atom.name);
        bool add_it = false;
        for (const auto& keep : keep_only_these_atoms) {
            if (atom_name == keep) {
                add_it = true;
                break;
            }
        }
        if (add_it) {
            addatom(atom_name,
                    atom.element.name(),
                    atom.pos.x, atom.pos.y, atom.pos.z,
                    std::string(1, atom.altloc),
                    atom.occ, atom.b_iso);
        }
    }
}

inline residue::residue(const gemmi::Residue& residue_ref,
                        const std::vector<std::string>& keep_only_these_atoms) {
    seqnum = residue_ref.seqid.num.value;
    ins_code = normalize_ins_code(residue_ref.seqid.icode);
    name = residue_ref.name;
    
    for (const auto& atom : residue_ref.atoms) {
        // Pad atom name to 4 characters in PDB format (1 leading space + name + trailing spaces)
        std::string atom_name = pad_atom_name_pdb(atom.name);
        bool add_it = false;
        for (const auto& keep : keep_only_these_atoms) {
            if (atom_name == keep) {
                add_it = true;
                break;
            }
        }
        if (add_it) {
            addatom(atom_name,
                    atom.element.name(),
                    atom.pos.x, atom.pos.y, atom.pos.z,
                    std::string(1, atom.altloc),
                    atom.occ, atom.b_iso);
        }
    }
}

inline void residue::addatom(std::string atom_name, std::string element,
                             float x, float y, float z, const std::string& altloc, float bf, float occupancy) {
    atom a(atom_name, element, x, y, z, altloc, occupancy, bf);
    atoms.push_back(a);
}

} // namespace minimol
} // namespace coot