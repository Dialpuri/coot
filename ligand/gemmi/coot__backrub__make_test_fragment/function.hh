#pragma once

#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <clipper/core/coords.h>

namespace coot {
namespace minimol {

struct atom {
    std::string name;
    clipper::Coord_orth pos;
    atom() : pos(0, 0, 0) {}
    atom(const std::string& n, double x, double y, double z)
        : name(n), pos(x, y, z) {}
};

class residue {
public:
    int seqnum;
    std::string ins_code;
    std::string name;
    std::vector<atom> atoms;

    residue() : seqnum(0), ins_code(""), name("") {}
    explicit residue(int sn) : seqnum(sn), ins_code(""), name("") {}
    residue(int sn, const std::string& resname) : seqnum(sn), ins_code(""), name(resname) {}

    // Build from gemmi residue, filtering to only the given atom names (trimmed)
    static residue from_gemmi(const gemmi::Residue* res, const std::vector<std::string>& include_atoms);

    unsigned int n_atoms() const { return static_cast<unsigned int>(atoms.size()); }
    const atom& operator[](unsigned int i) const { return atoms[i]; }
    atom& operator[](unsigned int i) { return atoms[i]; }
    bool is_empty() const { return atoms.empty(); }
};

class fragment {
    int residues_offset;
public:
    std::string fragment_id;
    std::vector<residue> residues;

    fragment() : residues_offset(0), fragment_id("") {
    }
    explicit fragment(const std::string& frag_id_in) : residues_offset(0), fragment_id(frag_id_in) {
    }

    void addresidue(const residue& res, int offset) {
        int idx = res.seqnum + offset;
        if (idx < 0) idx = 0;
        if (idx >= (int)residues.size()) {
            residues.resize(idx + 1);
        }
        residues[idx] = res;
    }

    int min_res_no() const { return residues_offset + 1; }
    int max_residue_number() const {
        if (residues.empty()) return residues_offset;
        return residues_offset + (int)residues.size() - 1;
    }

    const residue& operator[](int i) const {
        return residues[i - residues_offset];
    }
    residue& operator[](int i) {
        return residues[i - residues_offset];
    }

    int n_filled_residues() const {
        int count = 0;
        for (const auto& r : residues) {
            if (!r.is_empty()) count++;
        }
        return count;
    }
};

} // namespace minimol

namespace util {

inline clipper::Coord_orth rotate_around_vector(
    const clipper::Coord_orth& axis,
    const clipper::Coord_orth& pt,
    const clipper::Coord_orth& center,
    double angle_rad)
{
    clipper::Coord_orth v = pt - center;
    double c = std::cos(angle_rad);
    double s = std::sin(angle_rad);
    double one_minus_c = 1.0 - c;

    double len = std::sqrt(axis.x()*axis.x() + axis.y()*axis.y() + axis.z()*axis.z());
    if (len < 1e-10) {
        return pt;
    }
    clipper::Coord_orth k(axis.x()/len, axis.y()/len, axis.z()/len);

    clipper::Coord_orth kv(
        k.y()*v.z() - k.z()*v.y(),
        k.z()*v.x() - k.x()*v.z(),
        k.x()*v.y() - k.y()*v.x()
    );
    double k_dot_v = k.x()*v.x() + k.y()*v.y() + k.z()*v.z();

    return clipper::Coord_orth(
        center.x() + v.x()*c + kv.x()*s + k.x()*k_dot_v*one_minus_c,
        center.y() + v.y()*c + kv.y()*s + k.y()*k_dot_v*one_minus_c,
        center.z() + v.z()*c + kv.z()*s + k.z()*k_dot_v*one_minus_c
    );
}

inline std::string int_to_string(int i) {
    return std::to_string(i);
}

} // namespace util
} // namespace coot

namespace coot {
namespace backrub {

coot::minimol::fragment
make_test_fragment_gemmi(
    const gemmi::Residue* this_residue,
    const gemmi::Residue* prev_residue,
    const gemmi::Residue* next_residue,
    const std::string& chain_id,
    clipper::Coord_orth ca_prev,
    clipper::Coord_orth ca_next,
    double rotation_angle);

} // namespace backrub
} // namespace coot