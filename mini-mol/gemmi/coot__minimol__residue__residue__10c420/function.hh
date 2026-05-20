#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/elem.hpp>

namespace coot {
namespace minimol {

struct atom {
    std::string name;
    std::string element;
    float x, y, z;
    std::string altLoc;
    float occupancy;
    float temperature_factor;
    
    atom(std::string atom_name, std::string ele, float x_, float y_, float z_, 
         const std::string & altloc_, float occupancy_, float dbf)
        : name(atom_name), element(ele), x(x_), y(y_), z(z_),
          altLoc(altloc_), occupancy(occupancy_), temperature_factor(dbf) {}
    
    atom() : x(0), y(0), z(0), occupancy(0), temperature_factor(0) {}
};

struct residue {
    int seqnum;
    char ins_code;
    std::string name;
    std::vector<atom> atoms;
    
    residue(gemmi::Residue* residue_p) {
        seqnum = residue_p->seqid.num.value;
        ins_code = residue_p->seqid.icode;
        name = residue_p->name;
        
        for (const auto& atom : residue_p->atoms) {
            if (atom.element != gemmi::Element("X")) { // skip dummy atoms
                // PDB atom names and elements are left-aligned within fixed-width columns
                // Atom name: columns 13-16 (4 chars), Element: columns 77-78 (2 chars)
                // gemmi stores unpadded strings, so we need to reconstruct the PDB format
                
                // Add leading space to atom name (PDB format: " N  ")
                std::string padded_name = atom.name;
                if (padded_name.length() < 3) {
                    padded_name = " " + padded_name;
                }
                // Pad to 4 characters
                while (padded_name.length() < 4) padded_name += " ";
                
                // Add leading space to element (PDB format: " N", " C", etc.)
                std::string base_elem = atom.element.name();
                std::string padded_elem = std::string(" ") + base_elem;
                // Pad to 2 characters
                while (padded_elem.length() < 2) padded_elem += " ";
                
                addatom(padded_name, padded_elem,
                        atom.pos.x, atom.pos.y, atom.pos.z,
                        std::string(1, atom.altloc),
                        atom.occ,
                        atom.b_iso);
            }
        }
    }
    
    void addatom(const std::string& atom_name, const std::string& ele,
                 float x, float y, float z,
                 const std::string& altloc, float occupancy, float dbf) {
        atoms.emplace_back(atom_name, ele, x, y, z, altloc, occupancy, dbf);
    }
    
    size_t n_atoms() const { return atoms.size(); }
    
    const atom& operator[](size_t i) const { return atoms[i]; }
    atom& operator[](size_t i) { return atoms[i]; }
};

} // namespace minimol
} // namespace coot