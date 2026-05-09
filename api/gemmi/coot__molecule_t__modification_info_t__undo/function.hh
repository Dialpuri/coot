#pragma once
#include <string>
#include <vector>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {

namespace modification_info_t {

struct save_info_t {
    gemmi::Structure mol;
    std::string backup_name;
    
    save_info_t() = default;
    save_info_t(const gemmi::Structure& m, const std::string& name) 
        : mol(m), backup_name(name) {}
    
    gemmi::Structure get_mol() const { return mol; }
};

struct modification_info_t {
    std::string mol_name;
    bool is_mmcif_flag;
    std::vector<save_info_t> save_info;
    int modification_index;
    
    modification_info_t() : mol_name(""), is_mmcif_flag(false), modification_index(0) {}
    
    // Returns message (empty on success)
    std::string make_backup(const gemmi::Structure& mol, const std::string& name) {
        save_info.push_back(save_info_t(mol, name));
        modification_index = static_cast<int>(save_info.size());
        return "";
    }
    
    // Returns a new gemmi::Structure (copy) or empty structure if no backup available
    // The gemmi version returns the backup structure directly (value type instead of pointer)
    gemmi::Structure undo_gemmi(const gemmi::Structure& current_mol) {
        int idx = modification_index - 1;
        
        // If at the end, make a backup first (to allow undoing to this state)
        if (modification_index == static_cast<int>(save_info.size())) {
            make_backup(current_mol, "undo");
        }
        
        modification_index = idx;
        if (modification_index < 0)
            modification_index = 0;
        
        if (idx >= 0 && idx < static_cast<int>(save_info.size())) {
            return save_info[idx].get_mol();
        }
        
        // No backup available
        return gemmi::Structure();
    }
};

} // namespace modification_info_t
} // namespace coot