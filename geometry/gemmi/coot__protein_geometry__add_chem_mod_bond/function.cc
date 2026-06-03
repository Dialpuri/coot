#include "function.hh"

namespace coot {
namespace protein_geometry {

void add_chem_mod_bond_gemmi(
    std::vector<std::vector<std::string>>& rows,
    const std::vector<std::string>& tags,
    std::map<std::string, chem_mod>& mods)
{
    auto find_col = [&](const std::string& tag) -> int {
        for (int i = 0; i < static_cast<int>(tags.size()); ++i) {
            if (tags[i] == tag) return i;
        }
        return -1;
    };

    auto get_str = [&](int row_idx, const std::string& tag) -> std::pair<std::string, int> {
        int col = find_col(tag);
        if (col < 0 || row_idx < 0 || row_idx >= static_cast<int>(rows.size()))
            return std::make_pair(std::string(), 1);
        if (static_cast<std::size_t>(col) >= rows[row_idx].size())
            return std::make_pair(std::string(), 1);
        std::string val = rows[row_idx][col];
        if (val.empty())
            return std::make_pair(std::string(), 1);
        return std::make_pair(val, 0);
    };

    for (int j = 0; j < static_cast<int>(rows.size()); j++) {
        int ierr_tot = 0;

        std::string mod_id;
        std::string function;
        std::string atom_id_1;
        std::string atom_id_2;
        std::string new_type;
        double new_value_dist = 0.0;
        double new_value_dist_esd = 0.0;

        auto p = get_str(j, "mod_id");
        ierr_tot += p.second; mod_id = p.first;

        p = get_str(j, "function");
        ierr_tot += p.second; function = p.first;

        p = get_str(j, "atom_id_1");
        ierr_tot += p.second; atom_id_1 = p.first;

        p = get_str(j, "atom_id_2");
        ierr_tot += p.second; atom_id_2 = p.first;

        p = get_str(j, "new_type");
        ierr_tot += p.second; new_type = p.first;

        // new_value_dist
        {
            int col = find_col("new_value_dist");
            if (col < 0) {
                ierr_tot += 1;
            } else {
                std::string sv = rows[j][col];
                if (sv.empty() || sv == ".") {
                    ierr_tot += 1;
                } else {
                    new_value_dist = std::stod(sv);
                }
            }
        }

        // new_value_dist_esd
        {
            int col = find_col("new_value_dist_esd");
            if (col < 0) {
                ierr_tot += 1;
            } else {
                std::string sv = rows[j][col];
                if (sv.empty() || sv == ".") {
                    ierr_tot += 1;
                } else {
                    new_value_dist_esd = std::stod(sv);
                }
            }
        }

        if (ierr_tot == 0) {
            chem_mod_bond cmb(
                function,
                atom_id_mmdb_expand(atom_id_1),
                atom_id_mmdb_expand(atom_id_2),
                new_type,
                new_value_dist,
                new_value_dist_esd);
            mods[mod_id].add_mod_bond(cmb);
        }
    }
}

} // namespace protein_geometry
} // namespace coot