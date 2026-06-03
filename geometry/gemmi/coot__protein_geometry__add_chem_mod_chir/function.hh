#pragma once

#include <string>
#include <vector>
#include <map>
#include <iostream>

// Lightweight mmCIF-loop-like structure for _chem_mod_chir rows.
struct ChemModChirLoop {
    std::vector<std::string> tags;
    std::vector<std::vector<std::string>> rows;

    int GetLoopLength() const { return static_cast<int>(rows.size()); }

    const char* GetString(const std::string& tag, int row, int& ierr) const {
        ierr = 0;
        int col = -1;
        for (int i = 0; i < static_cast<int>(tags.size()); ++i) {
            if (tags[i] == tag) { col = i; break; }
        }
        if (col < 0 || row < 0 || row >= static_cast<int>(rows.size())) {
            ierr = 1; return nullptr;
        }
        if (rows[row].size() <= static_cast<size_t>(col)) {
            ierr = 1; return nullptr;
        }
        return rows[row][col].c_str();
    }
};

// chem_mod_chir — stores one chiral modification specification.
namespace coot {

struct chem_mod_chir {
    std::string function;
    std::string atom_id_centre;
    std::string atom_id_1;
    std::string atom_id_2;
    std::string atom_id_3;
    int volume_sign;

    chem_mod_chir(const std::string& f,
                  const std::string& c,
                  const std::string& a1,
                  const std::string& a2,
                  const std::string& a3,
                  int vs)
        : function(f), atom_id_centre(c),
          atom_id_1(a1), atom_id_2(a2), atom_id_3(a3),
          volume_sign(vs) {}
};

struct chem_mod {
    std::vector<chem_mod_chir> chir_mods;
    void add_mod_chir(const chem_mod_chir& c) { chir_mods.push_back(c); }
};

// Pad an atom name to 4-char PDB column format.
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

inline int chiral_volume_string_to_chiral_sign(const std::string& s) {
    if (s == "positive") return 1;
    if (s == "negative") return -1;
    return 0;
}

void add_chem_mod_chir_gemmi(ChemModChirLoop& loop, std::map<std::string, coot::chem_mod>& mods) {
    for (int j = 0; j < loop.GetLoopLength(); j++) {
        int ierr_tot = 0;
        int ierr;

        std::string mod_id;
        std::string function;
        std::string atom_id_centre;
        std::string atom_id_1;
        std::string atom_id_2;
        std::string atom_id_3;
        std::string new_volume_sign;

        char* s;

        s = const_cast<char*>(loop.GetString("mod_id", j, ierr));
        ierr_tot += ierr;
        if (s) mod_id = s;

        s = const_cast<char*>(loop.GetString("function", j, ierr));
        ierr_tot += ierr;
        if (s) function = s;

        s = const_cast<char*>(loop.GetString("atom_id_centre", j, ierr));
        ierr_tot += ierr;
        if (s) atom_id_centre = s;

        s = const_cast<char*>(loop.GetString("atom_id_1", j, ierr));
        ierr_tot += ierr;
        if (s) atom_id_1 = s;

        s = const_cast<char*>(loop.GetString("atom_id_2", j, ierr));
        ierr_tot += ierr;
        if (s) atom_id_2 = s;

        s = const_cast<char*>(loop.GetString("atom_id_3", j, ierr));
        ierr_tot += ierr;
        if (s) atom_id_3 = s;

        s = const_cast<char*>(loop.GetString("new_volume_sign", j, ierr));
        ierr_tot += ierr;
        if (s) new_volume_sign = s;

        int volume_sign_int =
            coot::protein_geometry::chiral_volume_string_to_chiral_sign(new_volume_sign);

        if (ierr_tot == 0) {
            coot::chem_mod_chir cmc(function,
                coot::atom_id_mmdb_expand(atom_id_centre),
                coot::atom_id_mmdb_expand(atom_id_1),
                coot::atom_id_mmdb_expand(atom_id_2),
                coot::atom_id_mmdb_expand(atom_id_3),
                volume_sign_int);
            mods[mod_id].add_mod_chir(cmc);
        } else {
            std::cout << "oops in add_chem_mod_chir ierr_tot is "
                      << ierr_tot << std::endl;
        }
    }
}

} // namespace protein_geometry
} // namespace coot