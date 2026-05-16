#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <map>

namespace coot {

// Lightweight mmCIF-loop-like structure for _chem_mod_atom rows.
// Mirrors the original mmdb::mmcif::PLoop interface just enough for the port.
struct ChemModAtomLoop {
    std::vector<std::string> tags;          // column names (e.g. "mod_id", "function", …)
    std::vector<std::vector<std::string>> rows; // one vector per row

    int GetLoopLength() const { return (int)rows.size(); }

    const char* GetString(const std::string& tag, int row, int& ierr) const {
        ierr = 0;
        int col = -1;
        for (int i = 0; i < (int)tags.size(); ++i) {
            if (tags[i] == tag) { col = i; break; }
        }
        if (col < 0 || row < 0 || row >= (int)rows.size()) {
            ierr = 1; return nullptr;
        }
        if (rows[row].size() <= (size_t)col) {
            ierr = 1; return nullptr;
        }
        return rows[row][col].c_str();
    }

    int GetReal(double& out, const std::string& tag, int row) const {
        int ierr = 0;
        const char* s = GetString(tag, row, ierr);
        if (ierr || s == nullptr) return 1;
        std::string sv(s);
        if (sv.empty() || sv == ".") return 1;
        out = std::stod(sv);
        return 0;
    }
};

// chem_mod_atom — stores one modification-atoms specification.
struct chem_mod_atom {
    std::string function;
    std::string atom_id;
    std::string new_atom_id;
    std::string new_type_symbol;
    std::string new_type_energy;
    double new_partial_charge;

    chem_mod_atom(const std::string& f, const std::string& a,
                  const std::string& n, const std::string& t,
                  const std::string& e, double q)
        : function(f), atom_id(a), new_atom_id(n),
          new_type_symbol(t), new_type_energy(e), new_partial_charge(q) {}
};

// chem_mod — container for a set of chem_mod_atom entries keyed by mod_id.
struct chem_mod {
    std::vector<chem_mod_atom> mod_atoms;

    void add_mod_atom(const chem_mod_atom& cma) {
        mod_atoms.push_back(cma);
    }
};

// Pad an atom name to 4-char PDB column format (identical logic to coot::atom_id_mmdb_expand).
inline std::string atom_id_mmdb_expand(const std::string& atomname) {
    std::string r;
    int ilen = (int)atomname.length();

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

// Ported version: iterates over a gemmi-compatible _chem_mod_atom loop and
// populates the mods map in the same way as the original MMDB function.
void add_chem_mod_atom_gemmi(const ChemModAtomLoop& loop,
                             std::map<std::string, chem_mod>& mods);

} // namespace protein_geometry
} // namespace coot