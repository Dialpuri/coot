#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <memory>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/make_asc/gemmi/function.hh"

namespace coot {

// --- save_info_t (gemmi port of modification_info save type) ---
struct save_info_t {
    std::string file_name;
    std::string modification_info_string;
    std::unique_ptr<gemmi::Structure> structure;

    save_info_t() = default;
    save_info_t(const std::string& fn, const std::string& mis)
        : file_name(fn), modification_info_string(mis) {}

    const gemmi::Structure* get_mol_gemmi() const { return structure.get(); }
};

// --- atom_selection_container_t (gemmi-compatible) ---
class atom_selection_container_t {
public:
    gemmi::Structure selected_atoms;
    int n_selected_atoms;

    atom_selection_container_t() : n_selected_atoms(0) {}
    atom_selection_container_t(const atom_selection_container_t& other)
        : selected_atoms(other.selected_atoms), n_selected_atoms(other.n_selected_atoms) {}
    atom_selection_container_t& operator=(const atom_selection_container_t& other) {
        if (this != &other) {
            selected_atoms = other.selected_atoms;
            n_selected_atoms = other.n_selected_atoms;
        }
        return *this;
    }

    void clear_up() {
        selected_atoms = gemmi::Structure();
        n_selected_atoms = 0;
    }

    void set_from_structure(const gemmi::Structure& st) {
        selected_atoms = st;
        n_selected_atoms = 0;
        for (const auto& model : selected_atoms.models) {
            for (const auto& chain : model.chains) {
                for (const auto& res : chain.residues) {
                    n_selected_atoms += static_cast<int>(res.atoms.size());
                }
            }
        }
    }
};

// --- modification_info_t (gemmi-compatible) ---
struct modification_info_t {
    std::string mol_name;
    bool is_mmcif_flag;
    std::vector<save_info_t> save_info;
    int modification_index;

    modification_info_t() : mol_name(""), is_mmcif_flag(false), modification_index(0) {}

    void make_backup(const gemmi::Structure& mol, const std::string& name) {
        save_info_t info(name, "");
        info.structure = std::make_unique<gemmi::Structure>(mol);
        save_info.push_back(std::move(info));
        modification_index = static_cast<int>(save_info.size());
    }

    // Undo: if at end, save current state first (without advancing index),
    // then go back one slot
    const gemmi::Structure* undo(const gemmi::Structure& current_mol) {
        if (modification_index == static_cast<int>(save_info.size())) {
            save_info_t info("undo", "");
            info.structure = std::make_unique<gemmi::Structure>(current_mol);
            save_info.push_back(std::move(info));
            // Do NOT update modification_index here — it stays pointing to
            // the same "current" position, and save_info grew.
        }
        if (modification_index <= 0)
            return nullptr;
        modification_index--;
        return save_info[modification_index].get_mol_gemmi();
    }

    // Redo: advance one slot in history, return pointer or nullptr at end
    const gemmi::Structure* redo() {
        int idx = modification_index + 1;
        if (idx >= static_cast<int>(save_info.size()))
            return nullptr;
        if (idx >= 0) {
            const gemmi::Structure* result = save_info[idx].get_mol_gemmi();
            modification_index = idx;
            return result;
        }
        return nullptr;
    }
};

// --- molecule_t (gemmi-compatible) ---
class molecule_t {
public:
    gemmi::Structure* mol;
    modification_info_t modification_info;
    atom_selection_container_t atom_sel;

    molecule_t() : mol(nullptr) {}
};

// --- redo_gemmi: main function being ported ---
//! gemmi port of coot::molecule_t::redo()
//! Calls the underlying modification_info redo, rebuilds atom_sel on success,
//! and always returns 0 (status).
inline int redo_gemmi(molecule_t& mol_obj) {
    int status = 0;

    const gemmi::Structure* mol_new = mol_obj.modification_info.redo();

    if (!mol_new) {
        std::cout << "ERROR:: undo failed" << std::endl;
    } else {
        mol_obj.atom_sel.set_from_structure(*mol_new);
    }

    return status;
}

// --- undo_gemmi: companion function ---
inline int undo_gemmi(molecule_t& mol_obj) {
    int status = 0;
    (void)mol_obj.modification_info.undo(*mol_obj.mol);
    return status;
}

} // namespace coot
