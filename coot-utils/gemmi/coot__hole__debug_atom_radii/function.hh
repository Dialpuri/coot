#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "coot/geometry/protein-geometry.hh"
#include "coot/geometry/residue-and-atom-specs.hh"

namespace coot::hole {

struct vdw_radii_result {
    std::vector<std::vector<std::vector<double>>> radii;
};

inline std::string pad_atom_name(const std::string& name) {
    std::string padded = name;
    while (padded.size() < 4) {
        padded = padded + ' ';
    }
    return padded;
}

inline vdw_radii_result assign_vdw_radii_gemmi(
    gemmi::Model& model,
    const coot::protein_geometry& geom)
{
    int imol = 0;
    bool use_vdwH_flag = false;

    std::map<std::pair<std::string, std::string>, double> cached_radii;
    vdw_radii_result result;

    result.radii.reserve(model.chains.size());

    for (gemmi::Chain& chain : model.chains) {
        std::vector<std::vector<double>> chain_radii;
        chain_radii.reserve(chain.residues.size());

        for (gemmi::Residue& residue : chain.residues) {
            std::vector<double> res_radii;
            res_radii.reserve(residue.atoms.size());

            std::string residue_name = residue.name;

            for (gemmi::Atom& atom : residue.atoms) {
                std::string padded_atom_name = pad_atom_name(atom.name);

                std::pair<std::string, std::string> p(padded_atom_name, residue_name);
                auto cit = cached_radii.find(p);

                double radius;
                if (cit != cached_radii.end()) {
                    radius = cit->second;
                } else {
                    radius = geom.get_vdw_radius(padded_atom_name, residue_name, imol, use_vdwH_flag);
                    cached_radii[p] = radius;
                }

                if (radius > 0) {
                    res_radii.push_back(radius);
                } else {
                    double default_radius = 1.7;
                    std::string ele = atom.element.name();
                    if (ele == "N")
                        default_radius = 1.55;
                    if (ele == "O")
                        default_radius = 1.52;
                    if (ele == "H")
                        default_radius = 1.2;
                    res_radii.push_back(default_radius);
                }
            }
            chain_radii.push_back(std::move(res_radii));
        }
        result.radii.push_back(std::move(chain_radii));
    }

    return result;
}

inline void debug_atom_radii_gemmi(const gemmi::Model& model, const vdw_radii_result& radii) {
    for (size_t ic = 0; ic < model.chains.size(); ic++) {
        const auto& chain = model.chains[ic];
        for (size_t ir = 0; ir < chain.residues.size(); ir++) {
            const auto& residue = chain.residues[ir];
            for (size_t ia = 0; ia < residue.atoms.size(); ia++) {
                const auto& atom = residue.atoms[ia];

                // Build atom_spec_t from gemmi fields (no _gemmi variant exists)
                coot::atom_spec_t spec;
                spec.chain_id = chain.name;
                spec.res_no = residue.seqid.num.value;
                spec.ins_code = std::string(1, residue.seqid.icode);
                spec.atom_name = atom.name;
                spec.alt_conf = std::string(1, atom.altloc);

                double radius = 0.0;
                if (ic < radii.radii.size() && ir < radii.radii[ic].size() && ia < radii.radii[ic][ir].size()) {
                    radius = radii.radii[ic][ir][ia];
                }

                std::cout << "   " << spec << " with radius " << radius << std::endl;
            }
        }
    }
}

} // namespace coot::hole