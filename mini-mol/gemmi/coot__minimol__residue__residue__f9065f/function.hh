#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>

namespace coot {
namespace minimol {

// Gemmi wrapper for residue - mimics coot::minimol::residue interface
class residue {
public:
    // ---- _gemmi factory constructors (the ported "residue" constructor) ----

    // Default: seqnum set to mmdb::MinInt4 equivalent (-2147483647)
    // Public member variables (matching original coot::minimol::residue interface)
    int seqnum;
    std::string ins_code;
    std::string name;

    // ---- _gemmi factory constructors (the ported "residue" constructor) ----

    // Default: seqnum set to mmdb::MinInt4 equivalent (-2147483647)
    static residue residue_gemmi() {
        residue r;
        r.seqnum   = -2147483647;
        r.ins_code = "";
        r.name     = "";
        return r;
    }

    // Constructor with seqnum
    static residue residue_gemmi(int seqnum) {
        residue r;
        r.seqnum   = seqnum;
        r.ins_code = "";
        r.name     = "";
        return r;
    }

    // Constructor with seqnum and residue name
    static residue residue_gemmi(int seqnum, const std::string& resname) {
        residue r;
        r.seqnum   = seqnum;
        r.ins_code = "";
        r.name     = resname;
        return r;
    }

    // ---- Legacy constructors (kept for compatibility) ----

    explicit residue(int seqnum) : seqnum(seqnum), ins_code(""), name("") {}

    residue(int seqnum, const std::string& resname)
        : seqnum(seqnum), ins_code(""), name(resname) {}

    residue() : seqnum(-2147483647), ins_code(""), name("") {}

    // Access underlying gemmi residue
    gemmi::Residue& get_residue() { return res_; }
    const gemmi::Residue& get_residue() const { return res_; }

    std::vector<gemmi::Atom>& atoms() { return res_.atoms; }
    const std::vector<gemmi::Atom>& atoms() const { return res_.atoms; }

    int n_atoms() const { return static_cast<int>(res_.atoms.size()); }

    bool is_empty() const { return res_.atoms.empty(); }

    // ---- addatom ----

    void addatom(const std::string& atom_name, const std::string& element,
                 float x, float y, float z, const std::string& altloc,
                 float bf, float occupancy) {
        gemmi::Atom atom;
        atom.name = atom_name;
        std::string elem = element;
        if (!elem.empty() && elem[0] == ' ') {
            elem = elem.substr(1);
        }
        atom.element = gemmi::Element(elem);
        atom.pos = gemmi::Position(x, y, z);
        atom.occ = occupancy;
        atom.b_iso = bf;
        if (!altloc.empty() && altloc != " ") {
            atom.altloc = altloc[0];
        } else {
            atom.altloc = '\0';
        }
        res_.atoms.push_back(atom);
    }

    void addatom(const gemmi::Atom& atom) {
        res_.atoms.push_back(atom);
    }

private:
    gemmi::Residue res_;
};

} // namespace minimol
} // namespace coot
