#pragma once

#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>

// Port of coot::beam_in_linked_residue::lsq_fit to gemmi.
// Original MMDB version: bool lsq_fit(mmdb::Residue *ref_res, mmdb::Residue *matcher_res,
//                                     mmdb::Residue *mov_res, const std::vector<std::string> &,
//                                     const std::vector<std::string> &).
// We replace mmdb::Residue* with gemmi::Residue* and handle atom lookup via gemmi APIs.

// Helper to get atoms by name from a gemmi residue
inline std::vector<gemmi::Atom*> beam_in_linked_residue_get_atoms_gemmi(gemmi::Residue *res,
                                                  const std::vector<std::string> &names) {
    std::vector<gemmi::Atom*> atoms;
    for (const auto& name : names) {
        for (auto& atom : res->atoms) {
            if (atom.name == name) {
                atoms.push_back(&atom);
                break;
            }
        }
    }
    return atoms;
}

// lsq_fit gemmi port. Returns true on success, false on mismatched atom counts.
inline bool beam_in_linked_residue_lsq_fit_gemmi(gemmi::Residue *ref_res,
                                                 gemmi::Residue *matcher_res,
                                                 gemmi::Residue *mov_res,
                                                 const std::vector<std::string> &lsq_atom_names_ref,
                                                 const std::vector<std::string> &lsq_atom_names_match) {
    bool status = false;
    std::vector<gemmi::Atom*> va_1 = beam_in_linked_residue_get_atoms_gemmi(ref_res, lsq_atom_names_ref);
    std::vector<gemmi::Atom*> va_2 = beam_in_linked_residue_get_atoms_gemmi(matcher_res, lsq_atom_names_match);

    if (va_1.size() != lsq_atom_names_ref.size()) {
        std::cout << "Mismatch atoms length in gemmi port (reference) "
                  << va_1.size() << " need " << lsq_atom_names_ref.size()
                  << std::endl;
    } else {
        if (va_1.size() != va_2.size()) {
            std::cout << "Mismatch atoms length in gemmi port (matcher) "
                      << va_1.size() << " vs " << va_2.size() << std::endl;
        } else {
            int n = static_cast<int>(lsq_atom_names_ref.size());
            std::vector<clipper::Coord_orth> co_1(n);
            std::vector<clipper::Coord_orth> co_2(n);
            for (int iat = 0; iat < n; iat++) {
                co_1[iat] = clipper::Coord_orth(va_1[iat]->pos.x, va_1[iat]->pos.y, va_1[iat]->pos.z);
                co_2[iat] = clipper::Coord_orth(va_2[iat]->pos.x, va_2[iat]->pos.y, va_2[iat]->pos.z);
            }
            clipper::RTop_orth rtop(co_1, co_2);

            // Transform coordinates manually using the RTop_orth
            for (auto& atom : mov_res->atoms) {
                clipper::Coord_orth c(atom.pos.x, atom.pos.y, atom.pos.z);
                clipper::Coord_orth c_new = rtop * c;
                atom.pos.x = c_new.x();
                atom.pos.y = c_new.y();
                atom.pos.z = c_new.z();
            }

            status = true;
        }
    }
    return status;
}