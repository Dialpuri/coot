#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <set>

#include "geometry/protein-geometry.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__deep_copy_this_residue_add_chain/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__put_amino_acid_residue_atom_in_standard_order/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_contact_indices_from_restraints__a15507/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__chi_angles__change_by__051042/gemmi/function.hh"

namespace coot {

// Inline helper: coot::util::is_standard_amino_acid_name (no gemmi port exists)
static inline bool is_standard_amino_acid_name_gemmi(const std::string& resname) {
    static const std::set<std::string> standard_aa = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL"
    };
    return standard_aa.count(resname) > 0;
}

void do_180_degree_side_chain_flip_gemmi(
    gemmi::Residue* residue,
    const std::string& alt_conf,
    coot::protein_geometry* geom_p) {

    if (!residue) return;

    std::string resname = residue->name;

    int nth_chi = -1; // unset

    if (resname == "ARG") nth_chi = 5;
    if (resname == "ASP") nth_chi = 2;
    if (resname == "ASN") nth_chi = 2;
    if (resname == "CYS") nth_chi = 1;
    if (resname == "GLN") nth_chi = 3;
    if (resname == "GLU") nth_chi = 3;
    if (resname == "PHE") nth_chi = 2;
    if (resname == "HIS") nth_chi = 2;
    if (resname == "SER") nth_chi = 1;
    if (resname == "THR") nth_chi = 1;
    if (resname == "VAL") nth_chi = 1;
    if (resname == "TRP") nth_chi = 2;
    if (resname == "TYR") nth_chi = 2;

    if (nth_chi != -1) {
        // Deep copy this residue (creates atoms with alt_conf)
        auto [residue_copy, new_chain] = coot::util::deep_copy_this_residue_add_chain_gemmi(
            gemmi::CRA{nullptr, residue, nullptr}, alt_conf, false, false);

        if (!residue_copy) return;

        // Check if we need reordering for standard amino acids
        if (is_standard_amino_acid_name_gemmi(resname)) {
            bool needs_reordering = false;
            int idx_N  = -1;
            int idx_CA = -1;
            int idx_CB = -1;

            int iat = 0;
            for (const auto& atom : residue_copy->atoms) {
                std::string at_name(atom.name);
                if (at_name == " N  ") idx_N  = iat;
                if (at_name == " CA ") idx_CA = iat;
                if (at_name == " CB ") idx_CB = iat;
                iat++;
            }

            if (idx_N != -1 && idx_CA != -1 && idx_N > idx_CA) {
                needs_reordering = true;
            }
            if (idx_CB != -1 && idx_CA != -1 && idx_CA > idx_CB) {
                needs_reordering = true;
            }

            if (needs_reordering) {
                coot::put_amino_acid_residue_atom_in_standard_order_gemmi(*residue_copy);
            }
        }

        // Use the gemmi ported change_by to rotate chi angle by 180 degrees
        double diff = 180.0;
        std::pair<short int, float> istat = coot::chi_angles::change_by_gemmi(
            nth_chi, diff, residue_copy, geom_p);

        if (istat.first) { // failure
            std::cout << "Failure to flip" << std::endl;
            delete residue_copy;
            return;
        }

        // Transfer coordinates of altconf atoms from residue_copy to residue
        // Normalize: MMDB "" and gemmi '\0' both mean "no alt loc"
        // When alt_conf is empty, we match atoms with altloc '\0' (the default)
        int nResidueAtoms = static_cast<int>(residue->atoms.size());
        int n_atom_residue_copy = static_cast<int>(residue_copy->atoms.size());

        char alt_conf_char = alt_conf.empty() ? '\0' : alt_conf[0];

        for (int iatc = 0; iatc < n_atom_residue_copy; iatc++) {
            const auto& atom_copy = residue_copy->atoms[iatc];
            if (atom_copy.altloc == alt_conf_char ||
                (atom_copy.altloc == '\0' && alt_conf.empty())) {
                std::string atom_copy_name = atom_copy.name;
                for (int iato = 0; iato < nResidueAtoms; iato++) {
                    std::string orig_atom_name = residue->atoms[iato].name;
                    if (orig_atom_name == atom_copy_name) {
                        residue->atoms[iato].pos.x = atom_copy.pos.x;
                        residue->atoms[iato].pos.y = atom_copy.pos.y;
                        residue->atoms[iato].pos.z = atom_copy.pos.z;
                    }
                }
            }
        }

        // Clean up the copy
        delete residue_copy;
    }
}

} // namespace coot
