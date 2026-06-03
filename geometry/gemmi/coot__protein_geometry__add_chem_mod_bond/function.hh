#pragma once

#include <string>
#include <vector>
#include <map>

namespace coot {

enum CHEM_MOD_FUNCTION {
    CHEM_MOD_FUNCTION_ADD=0,
    CHEM_MOD_FUNCTION_DELETE=1,
    CHEM_MOD_FUNCTION_MODIFY=2
};

struct chem_mod_bond {
    int function;
    std::string atom_id_1;
    std::string atom_id_2;
    std::string new_type;
    double new_value_dist;
    double new_value_dist_esd;

    chem_mod_bond(
        const std::string& func_str,
        const std::string& a1,
        const std::string& a2,
        const std::string& nt,
        double dist,
        double dist_esd)
        : atom_id_1(a1), atom_id_2(a2), new_type(nt),
          new_value_dist(dist), new_value_dist_esd(dist_esd)
    {
        if (func_str == "add")
            function = CHEM_MOD_FUNCTION_ADD;
        else if (func_str == "delete")
            function = CHEM_MOD_FUNCTION_DELETE;
        else
            function = CHEM_MOD_FUNCTION_MODIFY;
    }
};

struct chem_mod {
    std::vector<chem_mod_bond> bond_mods;

    void add_mod_bond(const chem_mod_bond& cmb) {
        bond_mods.push_back(cmb);
    }
};

inline std::string atom_id_mmdb_expand(const std::string& atomname) {
    std::string r;
    int ilen = static_cast<int>(atomname.length());

    if (ilen == 4) return atomname;

    if (ilen == 1) {
        r = " ";
        r += atomname;
        r += "  ";
    } else {
        if (ilen == 2) {
            static const char* metals[] = {
                "MG","NA","LI","AL","SI","CL","SC","TI","CR","MN","FE","CO",
                "NI","CU","ZN","GA","AS","SE","BR","RB","SR","RE","OS","IR",
                "PT","AU","HG","PB","BI"
            };
            bool done = false;
            for (const char* m : metals) {
                if (atomname == m) { done = true; break; }
            }
            if (done) {
                r += atomname;
                r += "  ";
            } else {
                r = " ";
                r += atomname;
                r += " ";
            }
        } else {
            if (ilen == 3) {
                r = " ";
                r += atomname;
            } else {
                r = atomname;
            }
        }
    }
    return r;
}

namespace protein_geometry {

void add_chem_mod_bond_gemmi(
    std::vector<std::vector<std::string>>& rows,
    const std::vector<std::string>& tags,
    std::map<std::string, chem_mod>& mods);

} // namespace protein_geometry
} // namespace coot