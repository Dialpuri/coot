#pragma once

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "/lmb/home/jdialpuri/Development/coot-dev/coot/geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_residue/gemmi/function.hh"

namespace coot {

inline int replace_residue_gemmi(
    const std::string& residue_cid,
    const std::string& new_residue_type,
    int imol_enc,
    const coot::protein_geometry& geom,
    gemmi::Model& model)
{
    int status = 0;

    gemmi::Residue* residue_p = coot::cid_to_residue_gemmi(residue_cid, model);
    if (!residue_p) {
        return 0;
    }

    // Get restraints for the new residue type
    std::pair<bool, coot::dictionary_residue_restraints_t> rp_new =
        geom.get_monomer_restraints(new_residue_type, imol_enc);
    if (!rp_new.first) {
        return 0;
    }
    const auto& restraints_new = rp_new.second;

    // Get restraints for the current residue type
    std::string current_residue_type = residue_p->name;
    std::pair<bool, coot::dictionary_residue_restraints_t> rp_current =
        geom.get_monomer_restraints(current_residue_type, imol_enc);
    if (!rp_current.first) {
        return 0;
    }
    const auto& restraints_current = rp_current.second;

    // Re-implement mutate_by_overlap logic using gemmi types

    // Trim whitespace helper for matching
    auto trim = [](const std::string& s) {
        size_t start = s.find_first_not_of(" ");
        if (start == std::string::npos) return std::string("");
        size_t end   = s.find_last_not_of(" ");
        return s.substr(start, end - start + 1);
    };

    // Build set of atom names from the new residue type
    std::map<std::string, int> new_atom_names;
    for (const auto& da : restraints_new.atom_info) {
        new_atom_names[da.atom_id_4c] = 0;
    }

    // 1) Keep atoms from current residue that exist in new type
    std::vector<gemmi::Atom> kept_atoms;
    for (auto& atom : residue_p->atoms) {
        std::string trimmed_name = trim(atom.name);
        if (new_atom_names.find(trimmed_name) != new_atom_names.end()) {
            kept_atoms.push_back(atom);
            new_atom_names[trimmed_name] = 1;
        }
    }

    // 2) For atoms that exist in both, update position from restraints
    for (auto& atom : kept_atoms) {
        std::string trimmed_name = trim(atom.name);
        for (const auto& da : restraints_new.atom_info) {
            if (da.atom_id_4c == trimmed_name) {
                atom.name = da.atom_id_4c;
                if (da.pdbx_model_Cartn_ideal.first) {
                    atom.pos.x = da.pdbx_model_Cartn_ideal.second.x();
                    atom.pos.y = da.pdbx_model_Cartn_ideal.second.y();
                    atom.pos.z = da.pdbx_model_Cartn_ideal.second.z();
                }
                break;
            }
        }
    }

    // 3) Add atoms from new type that don't exist in current
    for (const auto& kv : new_atom_names) {
        if (kv.second == 0) {
            // This atom is in new type but not in current — create it
            for (const auto& da : restraints_new.atom_info) {
                if (da.atom_id_4c == kv.first) {
                    gemmi::Atom new_atom;
                    new_atom.name = da.atom_id_4c;

                    // Determine element from type_symbol
                    std::string elem_sym;
                    if (!da.type_symbol.empty()) {
                        elem_sym = std::string(1, da.type_symbol[0]);
                        new_atom.element = gemmi::find_element(elem_sym.c_str());
                    } else {
                        new_atom.element = gemmi::Element("?");
                    }

                    if (da.pdbx_model_Cartn_ideal.first) {
                        new_atom.pos = gemmi::Position(
                            da.pdbx_model_Cartn_ideal.second.x(),
                            da.pdbx_model_Cartn_ideal.second.y(),
                            da.pdbx_model_Cartn_ideal.second.z());
                    } else {
                        new_atom.pos = gemmi::Position(0, 0, 0);
                    }
                    new_atom.occ = 1.0;
                    new_atom.b_iso = 0.0;
                    kept_atoms.push_back(new_atom);
                    break;
                }
            }
        }
    }

    residue_p->atoms = std::move(kept_atoms);

    // Update residue name
    residue_p->name = new_residue_type;

    status = 1;
    return status;
}

} // namespace coot