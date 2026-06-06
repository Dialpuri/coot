#include "function.hh"
#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include <gemmi/model.hpp>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_residue/gemmi/function.hh"
#include "geometry/protein-geometry.hh"
#include "mini-mol/atom-quads.hh"
#include "mini-mol/mini-mol.hh"
#include "coot-utils/atom-tree.hh"

namespace {

coot::minimol::residue make_minimol_residue_from_gemmi(const gemmi::Residue& gemmi_res) {
    coot::minimol::residue min_res;
    min_res.name = gemmi_res.name;
    min_res.seqnum = gemmi_res.seqid.num.value;
    min_res.ins_code = gemmi_res.seqid.icode;
    for (const gemmi::Atom& gemmi_atom : gemmi_res.atoms) {
        // gemmi atom names are PDB-padded (e.g. " N  "); trim to match CIF restraints
        std::string atom_name_trimmed = gemmi_atom.name;
        size_t start = atom_name_trimmed.find_first_not_of(" \t");
        if (start != std::string::npos) {
            atom_name_trimmed = atom_name_trimmed.substr(start);
        }
        size_t end = atom_name_trimmed.find_last_not_of(" \t");
        if (end != std::string::npos) {
            atom_name_trimmed = atom_name_trimmed.substr(0, end + 1);
        }

        // Normalize altloc: null or space → empty string (matches alt_conf="")
        std::string altloc_str;
        if (gemmi_atom.altloc != '\0' && gemmi_atom.altloc != ' ') {
            altloc_str = std::string(1, gemmi_atom.altloc);
        }
        min_res.addatom(
            atom_name_trimmed,
            gemmi_atom.element.name(),
            static_cast<float>(gemmi_atom.pos.x),
            static_cast<float>(gemmi_atom.pos.y),
            static_cast<float>(gemmi_atom.pos.z),
            altloc_str,
            static_cast<float>(gemmi_atom.b_iso),
            static_cast<float>(gemmi_atom.occ)
        );
    }
    return min_res;
}

} // anonymous namespace

int coot::rotate_around_bond_gemmi(
    const std::string& residue_cid,
    const std::string& alt_conf,
    coot::atom_name_quad quad,
    double torsion_angle,
    coot::protein_geometry& geom,
    gemmi::Model& model) {

    int status = 0;
    double r = -999.9;

    gemmi::Residue* residue_p = coot::cid_to_residue_gemmi(residue_cid, model);
    if (residue_p) {
        std::string res_name(residue_p->name);
        std::pair<bool, coot::dictionary_residue_restraints_t> restraints_info =
            geom.get_monomer_restraints(res_name, 0);
        if (!restraints_info.first) {
            std::cout << "WARNING:: set_torsion: No restraints for " << res_name << std::endl;
        } else {
            try {
                coot::minimol::residue min_res = make_minimol_residue_from_gemmi(*residue_p);
                coot::atom_tree_t tree(restraints_info.second, min_res, alt_conf);

                r = tree.set_dihedral(quad.atom_name(0),
                                      quad.atom_name(1),
                                      quad.atom_name(2),
                                      quad.atom_name(3),
                                      torsion_angle);
            }
            catch (const std::runtime_error& rte) {
                std::cout << "in set_torsion:: set_dihedral() error: " << rte.what() << std::endl;
            }
        }
    } else {
        std::cout << "failed to find residue " << residue_cid << std::endl;
    }
    return status;
}