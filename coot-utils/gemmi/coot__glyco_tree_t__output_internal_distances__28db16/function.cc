#include "function.hh"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace coot {

static std::string pad_to_4(const std::string &name) {
    if (name.size() >= 4) return name.substr(0, 4);
    std::ostringstream oss;
    oss << std::setw(4) << std::left << name;
    return oss.str();
}

void output_internal_distances_gemmi(
    const gemmi::Residue &residue,
    const gemmi::Residue *parent_p,
    double dist_crit,
    std::ofstream &f,
    int model_number,
    const std::string &chain_id) {

    bool include_hydrogen_atoms = false;

    // self distances - iterate by index to get forward pairs only
    const auto &atoms = residue.atoms;
    for (size_t iat = 0; iat < atoms.size(); iat++) {
        const gemmi::Atom &atom_i = atoms[iat];
        if (!atom_i.element.is_hydrogen()) {
            if (include_hydrogen_atoms || !atom_i.element.is_hydrogen()) {
                clipper::Coord_orth pos_atom_i = co_gemmi(&atom_i);
                
                for (size_t jat = iat; jat < atoms.size(); jat++) {
                    if (iat != jat) {
                        const gemmi::Atom &atom_j = atoms[jat];
                        if (!atom_j.element.is_hydrogen()) {
                            if (include_hydrogen_atoms || !atom_j.element.is_hydrogen()) {
                                clipper::Coord_orth pos_atom_j = co_gemmi(&atom_j);
                                double d = clipper::Coord_orth::length(pos_atom_i, pos_atom_j);
                                if (d < dist_crit && d > 0) {
                                    std::string atom_name_i = pad_to_4(atom_i.padded_name());
                                    std::string atom_name_j = pad_to_4(atom_j.padded_name());
                                    
                                    coot::atom_spec_t spec_i(
                                        chain_id.c_str(),
                                        residue.seqid.num.value,
                                        residue.seqid.icode == ' ' ? "" : std::string(1, residue.seqid.icode).c_str(),
                                        atom_name_i.c_str(),
                                        atom_i.altloc == '\0' ? "" : std::string(1, atom_i.altloc).c_str());
                                    spec_i.model_number = model_number;
                                    
                                    coot::atom_spec_t spec_j(
                                        chain_id.c_str(),
                                        residue.seqid.num.value,
                                        residue.seqid.icode == ' ' ? "" : std::string(1, residue.seqid.icode).c_str(),
                                        atom_name_j.c_str(),
                                        atom_j.altloc == '\0' ? "" : std::string(1, atom_j.altloc).c_str());
                                    spec_j.model_number = model_number;
                                    
                                    f << " " << spec_i.format() << " " << spec_j.format() << " " << d << std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // self-parent distances
    if (parent_p) {
        for (size_t iat = 0; iat < atoms.size(); iat++) {
            const gemmi::Atom &atom_i = atoms[iat];
            if (!atom_i.element.is_hydrogen()) {
                if (include_hydrogen_atoms || !atom_i.element.is_hydrogen()) {
                    clipper::Coord_orth pos_atom_i = co_gemmi(&atom_i);
                    
                    for (size_t jat = 0; jat < parent_p->atoms.size(); jat++) {
                        const gemmi::Atom &atom_j = parent_p->atoms[jat];
                        if (!atom_j.element.is_hydrogen()) {
                            if (include_hydrogen_atoms || !atom_j.element.is_hydrogen()) {
                                clipper::Coord_orth pos_atom_j = co_gemmi(&atom_j);
                                double d = clipper::Coord_orth::length(pos_atom_i, pos_atom_j);
                                if (d < dist_crit && d > 0) {
                                    std::string atom_name_i = pad_to_4(atom_i.padded_name());
                                    std::string atom_name_j = pad_to_4(atom_j.padded_name());
                                    
                                    coot::atom_spec_t spec_i(
                                        chain_id.c_str(),
                                        residue.seqid.num.value,
                                        residue.seqid.icode == ' ' ? "" : std::string(1, residue.seqid.icode).c_str(),
                                        atom_name_i.c_str(),
                                        atom_i.altloc == '\0' ? "" : std::string(1, atom_i.altloc).c_str());
                                    spec_i.model_number = model_number;
                                    
                                    coot::atom_spec_t spec_j(
                                        chain_id.c_str(),
                                        parent_p->seqid.num.value,
                                        parent_p->seqid.icode == ' ' ? "" : std::string(1, parent_p->seqid.icode).c_str(),
                                        atom_name_j.c_str(),
                                        atom_j.altloc == '\0' ? "" : std::string(1, atom_j.altloc).c_str());
                                    spec_j.model_number = model_number;
                                    
                                    f << " " << spec_i.format() << " " << spec_j.format() << " " << d << std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

} // namespace coot