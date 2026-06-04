#pragma once

#include <vector>
#include <map>
#include <set>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <gemmi/contact.hpp>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__is_main_chain_p/gemmi/function.hh"

namespace coot {
namespace daca {

static bool is_standard_amino_acid_name(const std::string& name) {
    return name == "ALA" || name == "CYS" || name == "ASP" || name == "GLU" ||
           name == "PHE" || name == "GLY" || name == "HIS" || name == "ILE" ||
           name == "LYS" || name == "LEU" || name == "MET" || name == "ASN" ||
           name == "PRO" || name == "GLN" || name == "ARG" || name == "SER" ||
           name == "THR" || name == "VAL" || name == "TRP" || name == "TYR";
}

std::vector<std::pair<gemmi::CRA, float>>
solvent_exposure_old_version_v2_gemmi(gemmi::Structure* st, bool side_chain_only) {

    std::vector<std::pair<gemmi::CRA, float>> v;
    if (!st) return v;

    std::map<gemmi::Residue*, std::set<gemmi::Atom*>> residue_neighbouring_atoms;

    if (st->models.empty()) return v;

    gemmi::Model& model = st->models[0];

    // Build neighbor contacts using ContactSearch
    float max_dist = 5.7;

    gemmi::NeighborSearch ns(model, st->cell, max_dist);
    ns.populate(false); // exclude hydrogens

    gemmi::ContactSearch cs(max_dist);
    cs.ignore = gemmi::ContactSearch::Ignore::SameResidue;

    std::vector<gemmi::ContactSearch::Result> contacts = cs.find_contacts(ns);

    for (const auto& contact : contacts) {
        // Only consider asymmetric unit contacts (no symmetry mates)
        if (contact.image_idx != 0) continue;

        gemmi::Atom* atom1 = contact.partner1.atom;
        gemmi::Atom* atom2 = contact.partner2.atom;

        if (!atom1 || !atom2) continue;

        gemmi::Residue* res1 = contact.partner1.residue;
        gemmi::Residue* res2 = contact.partner2.residue;

        if (!res1 || !res2) continue;
        if (res1 == res2) continue;

        std::string res_name_1 = res1->name;
        std::string res_name_2 = res2->name;

        if (res_name_1 == "HOH") continue;
        if (res_name_2 == "HOH") continue;
        if (!is_standard_amino_acid_name(res_name_1)) continue;

        // Original: if not side_chain_only, always add; if side_chain_only, only add if atom1 is NOT main chain
        if (!side_chain_only) {
            residue_neighbouring_atoms[res1].insert(atom2);
        } else {
            if (!coot::is_main_chain_p_gemmi(*atom1, *res1)) {
                residue_neighbouring_atoms[res1].insert(atom2);
            }
        }

        // Handle reverse direction (original SeekContacts makes reverses too)
        if (is_standard_amino_acid_name(res_name_2)) {
            if (!side_chain_only) {
                residue_neighbouring_atoms[res2].insert(atom1);
            } else {
                if (!coot::is_main_chain_p_gemmi(*atom2, *res2)) {
                    residue_neighbouring_atoms[res2].insert(atom1);
                }
            }
        }
    }

    // Now collect results for all residues that have neighbors
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            auto it = residue_neighbouring_atoms.find(&res);
            if (it != residue_neighbouring_atoms.end()) {
                v.push_back({gemmi::CRA{&chain, &res, nullptr}, static_cast<float>(it->second.size())});
            }
        }
    }

    return v;
}

} // namespace daca
} // namespace coot