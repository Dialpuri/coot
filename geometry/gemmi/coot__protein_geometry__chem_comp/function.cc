#include "function.hh"

#include <iostream>
#include <string>

namespace coot {
namespace protein_geometry {

// Stub implementations for linking — these would normally be provided by
// the coot library, but we need standalone compilation for the gemmi port.
void delete_mon_lib(const std::string& id, int imol_enc) {
    // stub — in coot, this deletes an existing mon_lib entry
}

void mon_lib_add_chem_comp(std::string id,
                           int imol_enc,
                           std::string three_letter_code,
                           std::string name,
                           std::string group,
                           int n_atoms_all,
                           int n_atoms_nh,
                           std::string description_level,
                           bool overwrite) {
    // stub — in coot, this adds a chem_comp entry to the mon_lib
}

std::string chem_comp_gemmi(gemmi::cif::Loop& loop, int imol_enc) {
    std::string returned_chem_comp;

    // Column indices based on the order the columns are added
    // id=0, three_letter_code=1, name=2, group=3,
    // number_atoms_all=4, number_atoms_nh=5, release_status=6, desc_level=7
    int col_id = 0;
    int col_three_letter_code = 1;
    int col_name = 2;
    int col_group = 3;
    int col_number_atoms_all = 4;
    int col_number_atoms_nh = 5;
    int col_release_status = 6;
    int col_desc_level = 7;

    int num_rows = static_cast<int>(loop.length());

    for (int j = 0; j < num_rows; j++) {
        std::string id = loop.val(j, col_id);
        std::string three_letter_code = loop.val(j, col_three_letter_code);
        std::string name = loop.val(j, col_name);
        std::string group = loop.val(j, col_group);
        if (group == "L-PEPTIDE")
            group = "L-peptide";
        int number_atoms_all = 0;
        std::string v_all = loop.val(j, col_number_atoms_all);
        if (!v_all.empty()) number_atoms_all = std::stoi(v_all);
        int number_atoms_nh = 0;
        std::string v_nh = loop.val(j, col_number_atoms_nh);
        if (!v_nh.empty()) number_atoms_nh = std::stoi(v_nh);
        std::string description_level = "None";
        std::string v_desc = loop.val(j, col_desc_level);
        if (!v_desc.empty()) {
            description_level = v_desc;
        } else {
            description_level = ".";
        }

        delete_mon_lib(id, imol_enc);
        mon_lib_add_chem_comp(id, imol_enc,
                              three_letter_code, name,
                              group, number_atoms_all, number_atoms_nh,
                              description_level, true);
        returned_chem_comp = id;
    }
    return returned_chem_comp;
}

} // namespace protein_geometry
} // namespace coot