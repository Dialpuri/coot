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
    explicit residue(int seqnum) : ins_code(""), name("") {
        res_.seqid.num.value = seqnum;
    }
    
    residue(int seqnum, const std::string& resname) : ins_code(""), name(resname) {
        res_.seqid.num.value = seqnum;
    }
    
    residue() : ins_code(""), name("") {
        res_.seqid.num.value = -2147483648; // mmdb::MinInt4 equivalent
    }
    
    int seqnum() const { return res_.seqid.num.value; }
    std::string ins_code;
    std::string name;
    
    // Access underlying gemmi residue
    gemmi::Residue& get_residue() { return res_; }
    const gemmi::Residue& get_residue() const { return res_; }
    
    std::vector<gemmi::Atom>& atoms() { return res_.atoms; }
    const std::vector<gemmi::Atom>& atoms() const { return res_.atoms; }
    
    void addatom(const std::string& atom_name, const std::string& element,
                 float x, float y, float z, const std::string& altloc,
                 float bf, float occupancy) {
        gemmi::Atom atom;
        atom.name = atom_name;
        // coot uses " N" (with leading space) for element names, gemmi uses "N"
        // Strip leading space if present to match gemmi's format
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