#pragma once
#include <string>
#include <vector>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {

class molecule_t {
public:
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
        
        std::string make_backup(const gemmi::Structure& mol, const std::string& name) {
            save_info.push_back(save_info_t(mol, name));
            modification_index = static_cast<int>(save_info.size());
            return "";
        }
        
        gemmi::Structure undo_gemmi(const gemmi::Structure& current_mol) {
            int idx = modification_index - 1;
            
            if (modification_index == static_cast<int>(save_info.size())) {
                make_backup(current_mol, "undo");
            }
            
            modification_index = idx;
            if (modification_index < 0)
                modification_index = 0;
            
            if (idx >= 0 && idx < static_cast<int>(save_info.size())) {
                return save_info[idx].get_mol();
            }
            
            return gemmi::Structure();
        }
    };
    
    modification_info_t modification_info;
    bool really_do_backups;
    gemmi::Structure atom_sel;
    
    molecule_t() : really_do_backups(true) {}
    
    int undo() {
        int status = 0;
        gemmi::Structure mol_new = modification_info.undo_gemmi(atom_sel);
        if (mol_new.models.empty()) {
            std::cout << "ERROR:: undo failed" << std::endl;
            status = -1;
        } else {
            atom_sel = mol_new;
        }
        return status;
    }
};

} // namespace coot