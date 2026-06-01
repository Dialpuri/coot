#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <map>
#include "geometry/protein-geometry.hh"

namespace coot {
namespace daca {

/// Holds a pointer to a gemmi atom, its parent residue, and the atom type.
struct typed_atom_gemmi {
    gemmi::Atom* atom;
    gemmi::Residue* residue;
    std::string type;
};

/// Pad gemmi atom names to PDB 4-character format.
/// PDB atom names use the convention: one leading space, then the name,
/// then trailing spaces to fill the 4-char field.
///   "N"   -> " N  "
///   "CA"  -> " CA "
///   "CD1" -> " CD1"
///   "OXT" -> " OXT"
inline std::string pad_atom_name(const std::string& name) {
    if (name.size() >= 4) return name.substr(0, 4);
    return " " + name + std::string(4 - name.size() - 1, ' ');
}

/// Gemmi port of coot::daca::make_typed_atoms.
///
/// Walks every chain/residue/atom in *model_p*, builds a dictionary of
/// monomer restraints from *geom*, then returns a vector of (atom, residue,
/// type-energy) for every atom whose type can be resolved.
///
/// Special-cases atom name "N" (-> " N  ") → type "NH1" (mirrors the MMDB original).
inline std::vector<typed_atom_gemmi> make_typed_atoms_gemmi(
    gemmi::Model* model_p,
    const coot::protein_geometry& geom) {

    std::vector<typed_atom_gemmi> v;
    std::map<std::string, dictionary_residue_restraints_t> dictionary_map;

    if (model_p) {
        // ---- First pass: populate restraint dictionary per residue type ----
        for (gemmi::Chain& chain : model_p->chains) {
            for (gemmi::Residue& res : chain.residues) {
                std::string res_type = res.name;
                auto it = dictionary_map.find(res_type);
                if (it == dictionary_map.end()) {
                    auto restraints =
                        geom.get_monomer_restraints(res_type,
                                                    protein_geometry::IMOL_ENC_ANY);
                    if (restraints.first) {
                        dictionary_map[res_type] = restraints.second;
                    }
                }
            }
        }

        // ---- Second pass: build typed atom list ----
        for (gemmi::Chain& chain : model_p->chains) {
            for (gemmi::Residue& res : chain.residues) {
                if (!res.atoms.empty()) {
                    std::string res_type = res.name;
                    auto it = dictionary_map.find(res_type);
                    if (it != dictionary_map.end()) {
                        for (gemmi::Atom& atom : res.atoms) {
                            // gemmi gives trimmed names ("N"), but type_energy needs
                            // 4-char padded names (" N  ").
                            std::string atom_name_padded = pad_atom_name(atom.name);
                            const std::string type = it->second.type_energy(atom_name_padded);

                            // Special case for backbone N (mirrors MMDB original which checked " N  ")
                            if (atom_name_padded == " N  ") {
                                v.push_back({&atom, &res, "NH1"});
                            } else if (!type.empty()) {
                                v.push_back({&atom, &res, type});
                            }
                        }
                    }
                }
            }
        }
    }

    return v;
}

} // namespace daca
} // namespace coot