#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <gemmi/elem.hpp>
#include <string>
#include <vector>

namespace coot { namespace minimol {

struct atom {
    std::string name;
    std::string element;
    gemmi::Position pos;
    double temperature_factor;
    double occupancy;
    std::string altLoc;

    atom() : pos(0, 0, 0), temperature_factor(0), occupancy(0) {}
};

struct residue {
    std::string name;
    int seqnum;
    std::string ins_code;
    std::vector<coot::minimol::atom> atoms;

    residue() : seqnum(0), ins_code("") {}
    residue(int seq, const std::string &nm) : name(nm), seqnum(seq), ins_code("") {}

    void addatom(const std::string &atom_name,
                 const std::string &ele,
                 float x, float y, float z,
                 const std::string &altloc,
                 float dbf,
                 float occ) {
        coot::minimol::atom a;
        a.name   = atom_name;
        a.element = ele;
        a.pos = gemmi::Position(static_cast<double>(x),
                                static_cast<double>(y),
                                static_cast<double>(z));
        a.altLoc            = altloc;
        a.temperature_factor = dbf;
        a.occupancy         = occ;
        atoms.push_back(a);
    }
};

}} // namespace coot::minimol

namespace coot {

static gemmi::El element_name_to_el(const std::string &s) {
    if (s == "H") return gemmi::El::H;
    if (s == "C") return gemmi::El::C;
    if (s == "N") return gemmi::El::N;
    if (s == "O") return gemmi::El::O;
    if (s == "S") return gemmi::El::S;
    if (s == "P") return gemmi::El::P;
    if (s == "Ca") return gemmi::El::Ca;
    if (s == "Mg") return gemmi::El::Mg;
    if (s == "Zn") return gemmi::El::Zn;
    if (s == "Fe") return gemmi::El::Fe;
    if (s == "Na") return gemmi::El::Na;
    if (s == "Cl") return gemmi::El::Cl;
    if (s == "K")  return gemmi::El::K;
    if (s == "Se") return gemmi::El::Se;
    return gemmi::El::X;
}

inline gemmi::Residue GetResidue_gemmi(const coot::minimol::residue &res_in) {
    gemmi::Residue res;
    res.seqid = gemmi::SeqId(res_in.seqnum,
                             res_in.ins_code.empty() ? ' ' : static_cast<char>(res_in.ins_code[0]));
    res.name = res_in.name;

    for (const auto &mat : res_in.atoms) {
        gemmi::Atom atom;
        atom.name    = mat.name;
        atom.element = gemmi::Element(element_name_to_el(mat.element));
        atom.pos     = mat.pos;
        atom.occ     = static_cast<float>(mat.occupancy);
        atom.b_iso   = static_cast<float>(mat.temperature_factor);
        atom.altloc  = mat.altLoc.empty() ? '\0' : static_cast<char>(mat.altLoc[0]);
        res.atoms.push_back(atom);
    }
    return res;
}

} // namespace coot